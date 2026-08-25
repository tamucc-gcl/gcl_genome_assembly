#!/usr/bin/env python3
# ======================================================================================
# gfa_hap_coverage.py
# Repo location: py_scripts/gfa_hap_coverage.py
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
#   <label>.hap_coverage_check.tsv    reconstructed marginal h(i) vs panacus
# ======================================================================================

import argparse
import gzip
import os
import sys
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


# ---- main ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="per-haplotype coverage decomposition of a pangenome GFA")
    ap.add_argument("--gfa", required=True, help="GFA (plain or gzipped); the clip GFA is the right input")
    ap.add_argument("--label", required=True, help="output filename prefix (species / taxid)")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--min-node-len", type=int, default=0,
                    help="ignore nodes shorter than this (default 0 = keep all)")
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

    # 3. marginal cross-check against panacus hist -------------------------------------
    marg = M.sum(axis=0)
    with open(op(".hap_coverage_check.tsv"), "w") as out:
        out.write("# bp_reconstructed = sum over haplotypes / coverage_level; compare to panacus hist\n")
        out.write("coverage_level\tsum_over_haplotypes_bp\tbp_reconstructed\n")
        for c in range(1, nhap + 1):
            out.write("%d\t%d\t%d\n" % (c, int(marg[c]), int(marg[c] // c)))

    recon_total = int(sum(marg[c] // c for c in range(1, nhap + 1)))
    sys.stderr.write("[gfa_hap_coverage] private_bp=%d  reconstructed_graph_bp=%d  (S-line bp=%d)\n"
                     % (tot_private, recon_total, total_bp))


if __name__ == "__main__":
    main()
