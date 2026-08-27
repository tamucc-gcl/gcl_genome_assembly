#!/usr/bin/env python3
# ======================================================================================
# gfa_hap_coverage.py  (v3)
# Repo location: py_scripts/gfa_hap_coverage.py
#
# v3 adds a per-contig breakdown of the private column. The v1 code below it -- the two
# numpy passes that produced output matching panacus exactly at all 10 coverage levels on
# Spratelloides -- is UNCHANGED. by_contig() is a separate third pass with its own parsing
# rather than a refactor of the verified path, deliberately: the small duplication buys
# zero risk to a result that is already validated against an independent implementation.
#
# Decomposes a pangenome GFA's coverage histogram BY HAPLOTYPE. panacus `hist` gives
# h(i) = bp of graph sequence covered by exactly i haplotypes, but it is a marginal: it
# cannot say WHICH haplotype owns the private (i == 1) sequence. This script computes the
# full joint table
#
#     M[haplotype, i] = bp of graph sequence that this haplotype traverses
#                       and that is covered by exactly i haplotypes in total
#
# from which the private-sequence breakdown is column i == 1, and the marginal
# sum_h M[h, i] / i reproduces panacus's h(i) exactly (built-in cross-check).
#
# Method: one bitmask word-set per node (bit g <=> haplotype group g traverses the node),
# so a node is counted once per haplotype no matter how many times it is walked. Coverage
# is the popcount. Grouping is by the first two PanSN fields (sample, haplotype), which is
# exactly what `panacus --groupby-haplotype` does, so the numbers are comparable.
#
# Handles both GFA path flavours emitted by minigraph-cactus / vg convert:
#   P  <path_name>  1+,2-,3+        ...   (PanSN name: sample#hap#contig)
#   W  <sample> <hap> <seq> <s> <e> >1<2>3
#
# Usage:
#   gfa_hap_coverage.py --gfa graph.gfa --label <species> --outdir .
#
# Outputs (in outdir):
#   <label>.hap_coverage_matrix.tsv   haplotype, coverage_level, bp        (the joint table)
#   <label>.hap_private.tsv           per-haplotype private bp + shares    (the plot input)
#   <label>.hap_private_by_contig.tsv per-contig private bp + share of that contig
#   <label>.hap_coverage_check.tsv    reconstructed marginal h(i) vs panacus
# ======================================================================================

import argparse
import gzip
import os
import sys
from array import array
from collections import OrderedDict

try:
    import numpy as np
except ImportError:  # pragma: no cover
    sys.exit("gfa_hap_coverage.py requires numpy (conda-forge::numpy)")


# ---- io ------------------------------------------------------------------------------
def gopen(path):
    """Open plain or gzipped GFA in binary mode, decided by magic bytes not extension."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    return open(path, "rb", buffering=1 << 20)


# strip orientation / separator characters so a step list becomes whitespace-delimited ids
_TR_P = bytes.maketrans(b"+-,", b"   ")   # P line: "12+,13-,14+"
_TR_W = bytes.maketrans(b"><", b"  ")     # W line: ">12<13>14"


def hap_key_from_path_name(name):
    """PanSN sample#hap#contig -> b'sample#hap'. Names without '#' become their own group."""
    parts = name.split(b"#")
    if len(parts) >= 2:
        return parts[0] + b"#" + parts[1]
    return name


# ---- pass 1: node lengths + haplotype group set --------------------------------------
def first_pass(gfa):
    """
    One read of the GFA collecting (a) segment lengths and (b) the haplotype group set.
    Groups are enumerated up front so bit assignment is sorted and run-to-run stable.
    Returns (nlen indexed by node id, n_segments, non_integer_example, sorted groups).
    """
    ids, lens = [], []
    bad = None
    n = 0
    seen = OrderedDict()
    with gopen(gfa) as fh:
        for line in fh:
            k = line[:2]
            if k == b"S\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                name, seq = f[1], f[2]
                if seq == b"*":
                    ln = 0
                    for tag in f[3:]:
                        if tag.startswith(b"LN:i:"):
                            ln = int(tag[5:])
                            break
                else:
                    ln = len(seq)
                try:
                    nid = int(name)
                except ValueError:
                    if bad is None:
                        bad = name.decode("utf-8", "replace")
                    continue
                ids.append(nid)
                lens.append(ln)
                n += 1
            elif k == b"P\t":
                f = line.split(b"\t", 2)
                if len(f) >= 2:
                    seen.setdefault(hap_key_from_path_name(f[1]), None)
            elif k == b"W\t":
                f = line.split(b"\t", 3)
                if len(f) >= 3:
                    seen.setdefault(f[1] + b"#" + f[2], None)
    groups = sorted(seen.keys())
    if n == 0:
        return None, 0, bad, groups
    ids = np.fromiter(ids, dtype=np.int64, count=n)
    lens = np.fromiter(lens, dtype=np.int64, count=n)
    nlen = np.zeros(int(ids.max()) + 1, dtype=np.int64)
    nlen[ids] = lens
    return nlen, n, bad, groups


# ---- pass 2: haplotype -> node bitmask ----------------------------------------------
def step_ids(field, translator, max_node):
    """Parse a step field into an int64 array of node ids, dropping out-of-range ids."""
    toks = field.translate(translator).split()
    if not toks:
        return None
    ids = np.array(toks, dtype="S").astype(np.int64)
    if ids.size and (ids.max() > max_node or ids.min() < 0):
        ids = ids[(ids >= 0) & (ids <= max_node)]
    return ids if ids.size else None


def build_masks(gfa, n_nodes, groups):
    """
    Second pass: bit g of node n is set iff haplotype group g traverses n. Set semantics,
    so repeated traversals by the same haplotype are counted once.
    Returns mask array of shape (n_words, n_nodes), plus P/W line counts.
    """
    gidx = {g: i for i, g in enumerate(groups)}
    n_words = (len(groups) + 63) // 64
    mask = np.zeros((n_words, n_nodes), dtype=np.uint64)
    n_p = n_w = 0
    max_node = n_nodes - 1
    with gopen(gfa) as fh:
        for line in fh:
            k = line[:2]
            if k == b"P\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                if len(f) < 3:
                    continue
                g = gidx.get(hap_key_from_path_name(f[1]))
                ids = step_ids(f[2], _TR_P, max_node)
                n_p += 1
            elif k == b"W\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                if len(f) < 7:
                    continue
                g = gidx.get(f[1] + b"#" + f[2])
                ids = step_ids(f[6], _TR_W, max_node)
                n_w += 1
            else:
                continue
            if g is None or ids is None:
                continue
            # duplicate ids within a line are harmless: every write stores (old | bit)
            mask[g // 64][ids] |= np.uint64(1) << np.uint64(g % 64)
    return mask, n_p, n_w


def popcount(mask):
    """Per-node haplotype coverage from the bitmask words."""
    n_words, n_nodes = mask.shape
    cov = np.zeros(n_nodes, dtype=np.int32)
    for w in range(n_words):
        word = mask[w]
        for b in range(64):
            shift = np.uint64(b)
            cov += ((word >> shift) & np.uint64(1)).astype(np.int32)
    return cov


# ---- pass 3: per-contig attribution of the private column ------------------------------
# Self-contained: parses P/W lines independently rather than sharing helpers with the two
# passes above, so nothing on the verified path is touched.
def by_contig(gfa, cov, nlen, groups):
    """
    Each P/W line is one contig walk, so `last` (the contig index that most recently touched
    a node) both de-duplicates repeats within a contig and, for private nodes, detects a
    second contig of the same haplotype walking the same node.

    Returns (order, priv, tot, multi_bp). A private node is credited to the FIRST contig that
    claims it; multi_bp is how much private sequence a second contig also walks -- a
    diagnostic, NOT additional bp, since it is already counted once in priv.

    tot can exceed a haplotype's hap_bp when two of its contigs share nodes; that is correct
    for "bp this contig walks", and is why pct_of_contig divides by tot.
    """
    # numpy scalar indexing is slow in a 1e8-iteration Python loop; copy to array module
    cov_l = cov.tolist()                      # values 0..nhap, all interned small ints
    nlen_a = array("q")
    if nlen_a.itemsize != 8:
        sys.exit("array('q') is not 8 bytes on this platform")
    nlen_a.frombytes(nlen.astype(np.int64).tobytes())

    cidx, order, priv, tot = {}, [], [], []
    last = array("i", bytes(4 * len(cov_l)))
    for i in range(len(last)):
        last[i] = -1
    multi_bp = 0
    known = set(groups)
    top = len(cov_l) - 1
    with gopen(gfa) as fh:
        for line in fh:
            k = line[:2]
            if k == b"P\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                if len(f) < 3:
                    continue
                p = f[1].split(b"#")
                if len(p) >= 3:
                    hap, contig = p[0] + b"#" + p[1], b"#".join(p[2:])
                elif len(p) == 2:
                    hap, contig = p[0] + b"#" + p[1], b"."
                else:
                    hap, contig = f[1], b"."
                fld, tr = f[2], _TR_P
            elif k == b"W\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                if len(f) < 7:
                    continue
                hap, contig, fld, tr = f[1] + b"#" + f[2], f[3], f[6], _TR_W
            else:
                continue
            if hap not in known:
                continue
            key = (hap, contig)
            c = cidx.get(key)
            if c is None:
                c = len(order)
                cidx[key] = c
                order.append(key)
                priv.append(0)
                tot.append(0)
            for i in map(int, fld.translate(tr).split()):
                if not (0 <= i <= top):
                    continue
                prev = last[i]
                if prev == c:
                    continue
                last[i] = c
                L = nlen_a[i]
                tot[c] += L
                if cov_l[i] == 1:
                    if prev == -1:
                        priv[c] += L
                    else:
                        multi_bp += L
    return order, priv, tot, multi_bp


# ---- pass 5: private segments (additive; passes 1-4 above are unchanged) ------------
def private_segments(gfa, cov, nlen, groups):
    """Maximal runs of consecutive private (cov == 1) nodes along each haplotype walk.

    Modelled on by_contig, which is the verified template for iterating P/W step lists, but
    it accumulates a POSITION along the walk instead of a bp total, so runs can be emitted
    with coordinates.

    Yields one record per run: (haplotype, contig, seg_index, n_nodes, seg_bp, start, end)
    where start/end are offsets in the contig's own coordinate space.

    Two deliberate differences from by_contig:

      1. No `last` deduplication. A private node walked twice is two positions in the
         contig, and a size spectrum is about positions. The bp a haplotype re-walks is
         accumulated separately as repeat_traversed_bp so the discrepancy against
         by_contig's deduplicated total is visible rather than mysterious.

      2. Offsets. W lines carry the contig interval in fields 5-6, so a clip-graph subpath
         starts where it actually starts. P lines have no coordinate and start at 0.
    """
    cov_l = cov.tolist()
    nlen_a = array("q")
    if nlen_a.itemsize != 8:
        sys.exit("array('q') is not 8 bytes on this platform")
    nlen_a.frombytes(nlen.astype(np.int64).tobytes())

    segs = []
    known = set(groups)
    top = len(cov_l) - 1
    # Marker array, not a per-haplotype reset: store WHICH haplotype last walked each
    # private node, exactly as by_contig stores which contig did. Zeroing a 153M-entry
    # array once per haplotype would be 10 x 153M Python iterations.
    seen_hap = array("i", bytes(4 * len(cov_l)))
    for _i in range(len(seen_hap)):
        seen_hap[_i] = -1
    hap_ids = {}
    repeat_bp = 0
    n_oor = 0

    with gopen(gfa) as fh:
        for line in fh:
            k = line[:2]
            if k == b"P\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                if len(f) < 3:
                    continue
                p = f[1].split(b"#")
                if len(p) >= 3:
                    hap, contig = p[0] + b"#" + p[1], b"#".join(p[2:])
                elif len(p) == 2:
                    hap, contig = p[0] + b"#" + p[1], b"."
                else:
                    hap, contig = f[1], b"."
                fld, tr, base = f[2], _TR_P, 0
            elif k == b"W\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                if len(f) < 7:
                    continue
                hap, contig, fld, tr = f[1] + b"#" + f[2], f[3], f[6], _TR_W
                try:
                    base = int(f[4])
                except ValueError:
                    base = 0
            else:
                continue
            if hap not in known:
                continue

            h = hap_ids.get(hap)
            if h is None:
                h = len(hap_ids)
                hap_ids[hap] = h

            hs = hap.decode("utf-8", "replace")
            cs = contig.decode("utf-8", "replace")
            pos = base
            run_bp = run_n = 0
            run_start = base
            idx = 0
            for i in map(int, fld.translate(tr).split()):
                if not (0 <= i <= top):
                    n_oor += 1
                    if run_n:
                        segs.append((hs, cs, idx, run_n, run_bp, run_start, run_start + run_bp))
                        idx += 1
                        run_bp = run_n = 0
                    continue
                L = nlen_a[i]
                if cov_l[i] == 1:
                    if seen_hap[i] == h:
                        repeat_bp += L
                    else:
                        seen_hap[i] = h
                    if run_n == 0:
                        run_start = pos
                    run_bp += L
                    run_n += 1
                elif run_n:
                    segs.append((hs, cs, idx, run_n, run_bp, run_start, run_start + run_bp))
                    idx += 1
                    run_bp = run_n = 0
                pos += L
            if run_n:
                segs.append((hs, cs, idx, run_n, run_bp, run_start, run_start + run_bp))
    return segs, repeat_bp, n_oor


SEG_BINS = [0, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, float("inf")]


def seg_bin(x):
    for lo, hi in zip(SEG_BINS, SEG_BINS[1:]):
        if lo <= x < hi:
            return (">=%d" % lo) if hi == float("inf") else ("%d-%d" % (lo, hi))
    return "unbinned"


# ---- main ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="per-haplotype coverage decomposition of a pangenome GFA")
    ap.add_argument("--gfa", required=True, help="GFA (plain or gzipped); the clip GFA is the right input")
    ap.add_argument("--label", required=True, help="output filename prefix (species / taxid)")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--min-contig-private-bp", type=int, default=0,
                    help="omit contigs below this private bp from the per-contig table (default 0 = all)")
    ap.add_argument("--min-node-len", type=int, default=0,
                    help="ignore nodes shorter than this (default 0 = keep all)")
    ap.add_argument("--min-private-bp", type=int, default=1000,
                    help="size floor for the private-segment TABLE and BED (default 1000). "
                         "The spectrum and totals are reported unfiltered as well, so the "
                         "floor never hides sequence -- it only bounds the row count.")
    ap.add_argument("--no-private-segments", action="store_true",
                    help="skip pass 5; the four original passes are unaffected either way")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    op = lambda suf: os.path.join(a.outdir, a.label + suf)

    nlen, n_seg, bad, groups = first_pass(a.gfa)
    if nlen is None:
        sys.exit("no S lines with integer names found in %s" % a.gfa)
    if bad is not None:
        sys.stderr.write("WARNING: non-integer segment name(s) skipped, e.g. %r\n" % bad)
    n_nodes = nlen.size
    total_bp = int(nlen.sum())
    sys.stderr.write("[gfa_hap_coverage] segments=%d  max_node_id=%d  graph_bp=%d\n"
                     % (n_seg, n_nodes - 1, total_bp))

    if not groups:
        sys.exit("no P or W lines found in %s" % a.gfa)
    names = [g.decode("utf-8", "replace") for g in groups]
    nhap = len(groups)
    sys.stderr.write("[gfa_hap_coverage] haplotype groups (n=%d): %s\n" % (nhap, ", ".join(names)))

    mask, n_p, n_w = build_masks(a.gfa, n_nodes, groups)
    sys.stderr.write("[gfa_hap_coverage] parsed %d P lines, %d W lines\n" % (n_p, n_w))

    cov = popcount(mask)
    keep = nlen >= max(a.min_node_len, 1)          # id 0 and absent ids have length 0
    n_orphan = int(np.count_nonzero(keep & (cov == 0)))
    if n_orphan:
        sys.stderr.write("WARNING: %d nodes carry sequence but no path (bp=%d)\n"
                         % (n_orphan, int(nlen[keep & (cov == 0)].sum())))

    # joint table M[hap, coverage_level] -------------------------------------------------
    M = np.zeros((nhap, nhap + 1), dtype=np.int64)
    for g in range(nhap):
        word = mask[g // 64]
        sel = (((word >> np.uint64(g % 64)) & np.uint64(1)) == 1) & keep
        if not sel.any():
            continue
        M[g] = np.bincount(cov[sel], weights=nlen[sel], minlength=nhap + 1)[: nhap + 1].astype(np.int64)

    hap_bp = M.sum(axis=1)
    private_bp = M[:, 1]
    tot_private = int(private_bp.sum())

    # 1. joint matrix ------------------------------------------------------------------
    with open(op(".hap_coverage_matrix.tsv"), "w") as out:
        out.write("# gfa=%s\n# haplotypes=%d\n# graph_bp=%d\n" % (os.path.basename(a.gfa), nhap, total_bp))
        out.write("haplotype\tcoverage_level\tbp\n")
        for g in range(nhap):
            for c in range(1, nhap + 1):
                out.write("%s\t%d\t%d\n" % (names[g], c, int(M[g, c])))

    # 2. private breakdown (the plot input) --------------------------------------------
    with open(op(".hap_private.tsv"), "w") as out:
        out.write("# total_private_bp=%d  n_haplotypes=%d  even_share_pct=%.4f\n"
                  % (tot_private, nhap, 100.0 / nhap))
        out.write("haplotype\tsample\thap\tprivate_bp\thap_bp\tpct_of_private\tpct_of_haplotype\n")
        order = np.argsort(-private_bp)
        for g in order:
            nm = names[g]
            smp, _, hp = nm.partition("#")
            pp = 100.0 * private_bp[g] / tot_private if tot_private else float("nan")
            ph = 100.0 * private_bp[g] / hap_bp[g] if hap_bp[g] else float("nan")
            out.write("%s\t%s\t%s\t%d\t%d\t%.6f\t%.6f\n"
                      % (nm, smp, hp, int(private_bp[g]), int(hap_bp[g]), pp, ph))

    # 3. per-contig attribution of the private column -----------------------------------
    order, cpriv, ctot, multi_bp = by_contig(a.gfa, cov, nlen, groups)
    hap_priv_lookup = dict(zip(groups, (int(x) for x in private_bp)))
    pct = lambda num, den: (100.0 * num / den) if den else None
    f4 = lambda v: "NA" if v is None else "%.4f" % v
    rows = []
    for c, (hg, cn) in enumerate(order):
        if cpriv[c] < a.min_contig_private_bp:
            continue
        rows.append((hg.decode("utf-8", "replace"), cn.decode("utf-8", "replace"),
                     cpriv[c], ctot[c], pct(cpriv[c], ctot[c]),
                     pct(cpriv[c], hap_priv_lookup.get(hg, 0)), pct(cpriv[c], tot_private)))
    rows.sort(key=lambda r: (r[0], -r[2]))
    with open(op(".hap_private_by_contig.tsv"), "w") as out:
        out.write("# private bp credited to the FIRST contig that walks each private node.\n")
        out.write("# multi_contig_private_bp=%d -- private sequence a SECOND contig of the same\n" % multi_bp)
        out.write("#   haplotype also walks. Diagnostic only: already counted once in private_bp.\n")
        out.write("# contig_graph_bp is bp this contig walks; contigs of one haplotype sharing a\n")
        out.write("#   node are each credited, so it can sum to more than that haplotype's hap_bp.\n")
        out.write("haplotype\tcontig\tprivate_bp\tcontig_graph_bp\tpct_of_contig\tpct_of_hap_private\tpct_of_pangenome_private\n")
        for r in rows:
            out.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\n" % (r[0], r[1], r[2], r[3], f4(r[4]), f4(r[5]), f4(r[6])))

    # 5. private segments: the size spectrum, and coordinates for everything downstream --
    if not a.no_private_segments:
        segs, repeat_bp, n_oor = private_segments(a.gfa, cov, nlen, groups)
        if n_oor:
            sys.stderr.write("WARNING: %d out-of-range node ids in step lists; each broke a "
                             "private run\n" % n_oor)

        seg_total = sum(s[4] for s in segs)
        kept = [s for s in segs if s[4] >= a.min_private_bp]
        kept_total = sum(s[4] for s in kept)

        with open(op(".private_segments.tsv"), "w") as out:
            out.write("# maximal runs of consecutive private (coverage==1) nodes along each\n")
            out.write("#   haplotype walk. start/end are offsets in the contig's own frame;\n")
            out.write("#   W lines supply a real start, P lines begin at 0.\n")
            out.write("# rows are filtered at --min-private-bp=%d; unfiltered totals are in\n"
                      % a.min_private_bp)
            out.write("#   the header below and in .private_segment_spectrum.tsv.\n")
            out.write("# segments_all=%d segment_bp_all=%d\n" % (len(segs), seg_total))
            out.write("# segments_kept=%d segment_bp_kept=%d\n" % (len(kept), kept_total))
            out.write("# repeat_traversed_bp=%d -- private bp a haplotype walks MORE THAN\n"
                      % repeat_bp)
            out.write("#   ONCE. by_contig (pass 4) credits such a node once, this pass counts\n")
            out.write("#   each position, so segment_bp_all exceeds hap_private.tsv by about\n")
            out.write("#   this amount. Not an error; it is repeat content within a haplotype.\n")
            out.write("haplotype\tcontig\tseg_index\tn_nodes\tseg_bp\tstart\tend\n")
            for s in sorted(kept, key=lambda x: (x[0], x[1], x[5])):
                out.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\n" % s)

        with open(op(".private_segments.bed"), "w") as out:
            out.write('track name="private_segments" description="private sequence >=%dbp"\n'
                      % a.min_private_bp)
            for s in sorted(kept, key=lambda x: (x[1], x[5])):
                out.write("%s\t%d\t%d\t%s|seg%d\t%d\t.\n"
                          % (s[1], s[5], s[6], s[0], s[2], min(1000, s[4])))

        # the histogram Chris asked for: count AND bp per bin, on the SV bin edges so the
        # two spectra can be read on one axis
        spec_n = {}
        spec_bp = {}
        per_hap_n = {}
        per_hap_bp = {}
        for hs, cs, si, nn, bp, st, en in segs:
            b = seg_bin(bp)
            spec_n[b] = spec_n.get(b, 0) + 1
            spec_bp[b] = spec_bp.get(b, 0) + bp
            per_hap_n[(hs, b)] = per_hap_n.get((hs, b), 0) + 1
            per_hap_bp[(hs, b)] = per_hap_bp.get((hs, b), 0) + bp
        with open(op(".private_segment_spectrum.tsv"), "w") as out:
            out.write("# UNFILTERED. Count and bp per size bin, on the same edges as the SV\n")
            out.write("#   spectrum so the two are directly comparable -- which is the actual\n")
            out.write("#   question: is the private-haplotype spectrum shaped like the SV one?\n")
            out.write("scope\tsize_bin\tn_segments\tsegment_bp\n")
            for b in sorted(spec_n, key=lambda x: SEG_BINS.index(
                    int(x.split("-")[0].lstrip(">="))) if x.split("-")[0].lstrip(">=").isdigit()
                    else 99):
                out.write("ALL\t%s\t%d\t%d\n" % (b, spec_n[b], spec_bp[b]))
            for (hs, b), n in sorted(per_hap_n.items()):
                out.write("%s\t%s\t%d\t%d\n" % (hs, b, n, per_hap_bp[(hs, b)]))

        sys.stderr.write("[gfa_hap_coverage] private segments: %d all (%d bp), "
                         "%d >=%dbp (%d bp)\n"
                         % (len(segs), seg_total, len(kept), a.min_private_bp, kept_total))
        sys.stderr.write("[gfa_hap_coverage]   repeat_traversed_bp=%d; segment_bp_all - "
                         "private_bp = %d\n" % (repeat_bp, seg_total - tot_private))

    # 4. marginal cross-check against panacus hist -------------------------------------
    marg = M.sum(axis=0)
    with open(op(".hap_coverage_check.tsv"), "w") as out:
        out.write("# bp_reconstructed = sum over haplotypes / coverage_level; compare to panacus hist\n")
        out.write("coverage_level\tsum_over_haplotypes_bp\tbp_reconstructed\n")
        for c in range(1, nhap + 1):
            out.write("%d\t%d\t%d\n" % (c, int(marg[c]), int(marg[c] // c)))

    recon_total = int(sum(marg[c] // c for c in range(1, nhap + 1)))
    sys.stderr.write("[gfa_hap_coverage] private_bp=%d  reconstructed_graph_bp=%d  (S-line bp=%d)\n"
                     % (tot_private, recon_total, total_bp))
    sys.stderr.write("[gfa_hap_coverage] contigs=%d  multi-contig private bp=%d (%.3f%% of private)\n"
                     % (len(order), multi_bp, 100.0 * multi_bp / tot_private if tot_private else 0.0))
    for r in sorted(rows, key=lambda x: -x[2])[:10]:
        sys.stderr.write("    %-26s %-18s %8.1f Mb private  %5.1f%% of contig (%.1f Mb)\n"
                         % (r[0], r[1], r[2] / 1e6, r[4] or 0.0, r[3] / 1e6))


if __name__ == "__main__":
    main()
