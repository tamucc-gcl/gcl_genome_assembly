#!/usr/bin/env python3
# ======================================================================================
# extract_private_fasta.py
#
# Emits one FASTA per haplotype containing that haplotype's PRIVATE segments -- maximal runs
# of consecutive nodes covered by exactly one haplotype -- at or above a size floor.
#
# WHY THIS RECOMPUTES COVERAGE INSTEAD OF READING THE BED
# ------------------------------------------------------
# gfa_hap_coverage.py already emits private_segments.bed, but in PATH coordinates
# (contig/start/end), not node ids. Turning a path interval back into sequence needs the node
# walk, which needs coverage, which is most of that script. Passing node id lists between the
# two processes would mean emitting ~85M ids as text. Recomputing is one extra GFA read and a
# popcount -- the HAP_COVERAGE task runs in minutes -- and it keeps this script independent,
# so a change to the BED format cannot silently break extraction.
#
# The coverage definition is deliberately identical to gfa_hap_coverage.py's: a node is
# private when exactly one haplotype group walks it. If the two ever disagree, the segment
# count and bp printed at the end will not match private_segment_spectrum.tsv, which is the
# intended tripwire.
#
# WHAT THE OUTPUT IS FOR
# ----------------------
# Two questions, both live:
#
#   1. Is "private" sequence real, or graph collapse? A segment that aligns cleanly to three
#      other assemblies is not private -- either the graph failed to merge it or the
#      alignment failed. Only mapping it back against the other assemblies can tell.
#   2. Is the elevated private fraction on chr8 (23.6%) and chr9 (22.9%) -- against an 11%
#      floor across other chromosomes, consistent across all ten haplotypes and both graph
#      flavours -- real divergence or collapsed repeat?
#
# Headers carry haplotype, contig and offset so mapping and k-mer results stay joinable to
# the segment table without a coordinate lookup.
#
# FLAVOUR MATTERS MORE THAN IT LOOKS
# ----------------------------------
# On the clip graph this yields 163,418 segments / 675 Mb; on the full graph 170,311 / 1371 Mb.
# Clipping removes 698 Mb of private sequence -- 99.1% of everything it removes -- and the
# entire difference lands in the >=1 kb segments that carry 74% of the signal. Run both.
#
# MEMORY
# ------
# Peak is roughly: node lengths (int32 over max node id) + coverage counts (uint8) + per
# segment node id lists + the sequences of nodes inside kept segments. On the Spratelloides
# full graph (153M nodes, 1.37 Gb of private sequence >=1 kb) that is ~3 GB.
#
# Dependencies: python3 + numpy.
#
# USAGE
#   extract_private_fasta.py --gfa graph.gfa.gz --label 373251.full --outdir . \
#       [--min-bp 1000] [--min-node-len 0]
# ======================================================================================

import argparse
import gzip
import os
import sys
from array import array

try:
    import numpy as np
except ImportError:
    sys.exit("extract_private_fasta.py requires numpy (conda-forge::numpy)")

# translate the orientation characters out of P/W step lists, leaving whitespace-separated ids
_TR_P = bytes.maketrans(b"+-,", b"   ")
_TR_W = bytes.maketrans(b"><", b"  ")


def gopen(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rb") if gz else open(path, "rb")


def hap_key(name):
    """PanSN path name -> haplotype group key.

    Sde-CBau_104#1#chr10_1  -> Sde-CBau_104#1
    Sde-CMat_203_hap2#0#chr8 -> Sde-CMat_203_hap2#0
    Matches gfa_hap_coverage.py so the two agree on what a haplotype is.
    """
    p = name.split(b"#")
    return p[0] + b"#" + p[1] if len(p) >= 2 else name


# --------------------------------------------------------------------------------------
def pass_lengths(gfa):
    """S lines -> (int32 lengths by node id, max_id, n_segments, total_bp)."""
    cap = 1 << 24
    arr = np.zeros(cap, dtype=np.int32)
    max_id = n = 0
    with gopen(gfa) as fh:
        for line in fh:
            if line[:2] != b"S\t":
                continue
            f = line.split(b"\t", 3)
            try:
                nid = int(f[1])
            except ValueError:
                continue
            if nid >= cap:
                while nid >= cap:
                    cap <<= 1
                g = np.zeros(cap, dtype=np.int32)
                g[:len(arr)] = arr
                arr = g
            arr[nid] = len(f[2].rstrip())
            n += 1
            if nid > max_id:
                max_id = nid
    arr = arr[:max_id + 1]
    return arr, max_id, n, int(arr.sum())


def pass_coverage(gfa, max_id, min_node_len, lengths):
    """P/W lines -> per-node count of DISTINCT haplotype groups walking it, and the group list.

    Counting distinct groups rather than steps is the point: a haplotype that walks a node
    twice must not make it look shared.
    """
    groups = {}
    # one bitmask byte-pair per node is prohibitive at 153M nodes x 10 haplotypes, so instead
    # track "last group seen" plus a count, exactly as gfa_hap_coverage.py does
    last = array("i", bytes(4 * (max_id + 1)))
    for i in range(len(last)):
        last[i] = -1
    cov = array("B", bytes(max_id + 1))

    with gopen(gfa) as fh:
        for line in fh:
            k = line[:2]
            if k == b"P\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                if len(f) < 3:
                    continue
                hk, fld, tr = hap_key(f[1]), f[2], _TR_P
            elif k == b"W\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                if len(f) < 7:
                    continue
                hk, fld, tr = f[1] + b"#" + f[2], f[6], _TR_W
            else:
                continue
            g = groups.setdefault(hk, len(groups))
            for i in map(int, fld.translate(tr).split()):
                if 0 <= i <= max_id and last[i] != g:
                    last[i] = g
                    if cov[i] < 255:
                        cov[i] += 1
    return cov, groups


def pass_segments(gfa, cov, lengths, max_id, min_bp):
    """P/W lines -> private runs >= min_bp, as (hap, contig, start, bp, [node ids]).

    Node ORDER is retained so sequence can be assembled without a third path walk. Runs are
    positional: a private node walked twice is two pieces of that contig.
    """
    segs = []
    with gopen(gfa) as fh:
        for line in fh:
            k = line[:2]
            if k == b"P\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                if len(f) < 3:
                    continue
                p = f[1].split(b"#")
                hk = hap_key(f[1])
                contig = b"#".join(p[2:]) if len(p) >= 3 else b"."
                fld, tr, base = f[2], _TR_P, 0
            elif k == b"W\t":
                f = line.rstrip(b"\r\n").split(b"\t")
                if len(f) < 7:
                    continue
                hk, contig, fld, tr = f[1] + b"#" + f[2], f[3], f[6], _TR_W
                try:
                    base = int(f[4])
                except ValueError:
                    base = 0
            else:
                continue

            pos = base
            run_ids, run_bp, run_start = [], 0, base
            for i in map(int, fld.translate(tr).split()):
                if not (0 <= i <= max_id):
                    if run_bp >= min_bp:
                        segs.append((hk, contig, run_start, run_bp, run_ids))
                    run_ids, run_bp = [], 0
                    continue
                L = int(lengths[i])
                if cov[i] == 1:
                    if not run_ids:
                        run_start = pos
                    run_ids.append(i)
                    run_bp += L
                else:
                    if run_bp >= min_bp:
                        segs.append((hk, contig, run_start, run_bp, run_ids))
                    run_ids, run_bp = [], 0
                pos += L
            if run_bp >= min_bp:
                segs.append((hk, contig, run_start, run_bp, run_ids))
    return segs


def pass_sequences(gfa, wanted):
    """S lines -> {node id: sequence bytes} for the given id set only."""
    seqs = {}
    with gopen(gfa) as fh:
        for line in fh:
            if line[:2] != b"S\t":
                continue
            f = line.split(b"\t", 3)
            try:
                nid = int(f[1])
            except ValueError:
                continue
            if nid in wanted:
                seqs[nid] = f[2].rstrip()
    return seqs


# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--label", required=True, help="e.g. 373251.full")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-bp", type=int, default=1000,
                    help="size floor on a private run (default 1000). 99.0%% of private "
                         "segments are under 100 bp and carry 7%% of the sequence, so a floor "
                         "removes almost all rows and almost no signal.")
    ap.add_argument("--min-node-len", type=int, default=0)
    ap.add_argument("--wrap", type=int, default=0, help="FASTA line wrap (0 = single line)")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)

    sys.stderr.write("[private_fasta] pass 1/4: node lengths\n")
    lengths, max_id, n_seg, graph_bp = pass_lengths(a.gfa)
    sys.stderr.write("[private_fasta]   %d segments, max id %d, graph %d bp\n"
                     % (n_seg, max_id, graph_bp))

    sys.stderr.write("[private_fasta] pass 2/4: coverage\n")
    cov, groups = pass_coverage(a.gfa, max_id, a.min_node_len, lengths)
    sys.stderr.write("[private_fasta]   %d haplotype groups: %s\n"
                     % (len(groups), ", ".join(sorted(g.decode() for g in groups))))

    sys.stderr.write("[private_fasta] pass 3/4: private runs >= %d bp\n" % a.min_bp)
    segs = pass_segments(a.gfa, cov, lengths, max_id, a.min_bp)
    tot_bp = sum(s[3] for s in segs)
    sys.stderr.write("[private_fasta]   %d segments, %d bp\n" % (len(segs), tot_bp))
    if not segs:
        sys.stderr.write("[private_fasta] nothing to extract; writing empty manifest\n")

    wanted = set()
    for _, _, _, _, ids in segs:
        wanted.update(ids)
    sys.stderr.write("[private_fasta] pass 4/4: sequences for %d nodes\n" % len(wanted))
    seqs = pass_sequences(a.gfa, wanted)
    missing = len(wanted) - len(seqs)
    if missing:
        sys.stderr.write("WARNING: %d node ids had no S-line sequence\n" % missing)

    # ---- write one FASTA per haplotype -----------------------------------------------
    handles, counts, bps = {}, {}, {}
    manifest = os.path.join(a.outdir, "%s.private_fasta_manifest.tsv" % a.label)
    for hk, contig, start, bp, ids in segs:
        h = hk.decode("utf-8", "replace")
        safe = h.replace("#", "_").replace("/", "_")
        if h not in handles:
            handles[h] = open(os.path.join(a.outdir, "%s.%s.private.fa" % (a.label, safe)), "w")
            counts[h] = bps[h] = 0
        seq = b"".join(seqs.get(i, b"") for i in ids)
        counts[h] += 1
        bps[h] += len(seq)
        # header carries everything needed to join mapping and k-mer results back to the
        # segment table without a coordinate lookup
        hdr = "%s|%s|%d|%d|seg%d" % (h, contig.decode("utf-8", "replace"),
                                     start, start + bp, counts[h])
        fh = handles[h]
        fh.write(">%s\n" % hdr)
        s = seq.decode("ascii", "replace")
        if a.wrap and a.wrap > 0:
            for i in range(0, len(s), a.wrap):
                fh.write(s[i:i + a.wrap] + "\n")
        else:
            fh.write(s + "\n")
    for fh in handles.values():
        fh.close()

    with open(manifest, "w") as out:
        out.write("# one FASTA per haplotype. Headers are\n")
        out.write("#   <haplotype>|<contig>|<start>|<end>|seg<N>\n")
        out.write("# so PRIVATE_MAP and PRIVATE_KMER results join back without a lookup.\n")
        out.write("# segment_bp is the private run length; fasta_bp is the sequence written.\n")
        out.write("# They differ only if an S-line was missing for a node (see warnings).\n")
        out.write("label\thaplotype\tfasta\tn_segments\tfasta_bp\n")
        for h in sorted(handles):
            safe = h.replace("#", "_").replace("/", "_")
            out.write("%s\t%s\t%s.%s.private.fa\t%d\t%d\n"
                      % (a.label, h, a.label, safe, counts[h], bps[h]))

    with open(os.path.join(a.outdir, "%s.private_fasta_audit.tsv" % a.label), "w") as out:
        out.write("metric\tvalue\n")
        out.write("label\t%s\n" % a.label)
        out.write("gfa\t%s\n" % os.path.basename(a.gfa))
        out.write("graph_segments\t%d\ngraph_bp\t%d\n" % (n_seg, graph_bp))
        out.write("haplotype_groups\t%d\n" % len(groups))
        out.write("min_bp\t%d\n" % a.min_bp)
        out.write("segments_kept\t%d\nsegment_bp\t%d\n" % (len(segs), tot_bp))
        out.write("nodes_needed\t%d\nnodes_missing_sequence\t%d\n" % (len(wanted), missing))
        out.write("fasta_bp_total\t%d\n" % sum(bps.values()))
        out.write("# segment_bp and fasta_bp_total must agree. They also must match the\n")
        out.write("# >=min_bp rows of private_segment_spectrum.tsv from gfa_hap_coverage.py --\n")
        out.write("# the two compute coverage independently, so a mismatch means one of them\n")
        out.write("# is wrong about what 'private' is.\n")

    sys.stderr.write("[private_fasta] wrote %d haplotype FASTAs, %d bp total\n"
                     % (len(handles), sum(bps.values())))
    if tot_bp != sum(bps.values()):
        sys.stderr.write("WARNING: segment_bp %d != fasta_bp %d -- missing node sequences\n"
                         % (tot_bp, sum(bps.values())))


if __name__ == "__main__":
    main()
