#!/usr/bin/env python3
# ======================================================================================
# summarise_private_map.py
#
# Turns a PAF of one haplotype's private segments aligned against all assemblies into a
# per-segment table answering: does this "private" sequence exist elsewhere?
#
# WHY THIS EXISTS
# ---------------
# "Private" is defined by GRAPH COVERAGE -- a node walked by exactly one haplotype. That is
# not the same as sequence absent from other assemblies. A segment can be private in the
# graph because the aligner failed to merge it, or because minigraph-cactus collapsed a
# repeat, while the sequence itself sits plainly in three other assemblies. Nothing else in
# the pipeline can tell those apart, and every private-sequence figure depends on the
# distinction.
#
# The live question it is aimed at: chr8 and chr9 carry 23.6% and 22.9% private sequence
# against an ~11% floor across other chromosomes, consistently across all ten haplotypes and
# both graph flavours. Real divergence, or collapsed repeat? If those segments align cleanly
# to other assemblies, it is collapse.
#
# TWO THINGS DONE CAREFULLY
# -------------------------
# 1. SELF-HITS ARE EXCLUDED BY ASSEMBLY TAG, NOT BY IDENTITY. A private segment is present in
#    its own assembly by construction and will align there at ~100%. Filtering on identity
#    would also discard genuine near-identical copies elsewhere -- exactly the collapsed-repeat
#    case we are trying to detect. Target names are <assembly_id>::<contig> and the query's
#    own assembly id is passed in explicitly, because it cannot be recovered from the PanSN
#    haplotype name (sampleOf() maps dots to underscores and is not invertible).
#
# 2. ALIGNED FRACTION USES MERGED QUERY INTERVALS. Summing PAF block lengths double-counts
#    overlapping alignments, which is the same error as measuring an alignment envelope
#    instead of a merged footprint -- it has appeared three times in this project already, so
#    it is merged here per (segment, target assembly) and again across assemblies.
#
# Dependencies: python3 stdlib.
#
# USAGE
#   summarise_private_map.py --paf x.paf.gz --self-assembly Sde-CBau_104_hap1 \
#       --haplotype 'Sde-CBau_104#1' --label 373251.full --outdir .
# ======================================================================================

import argparse
import gzip
import os
import sys
from collections import defaultdict


def popen(p):
    with open(p, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return (gzip.open if gz else open)(p, "rt", encoding="utf-8", errors="replace")


def merged_len(iv):
    """Union length of [start, end) intervals. NOT the sum of block lengths."""
    if not iv:
        return 0
    iv = sorted(iv)
    total, s, e = 0, iv[0][0], iv[0][1]
    for a, b in iv[1:]:
        if a > e:
            total += e - s
            s, e = a, b
        elif b > e:
            e = b
    return total + (e - s)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True)
    p.add_argument("--self-assembly", required=True,
                   help="assembly id of the query haplotype; its hits are EXCLUDED")
    p.add_argument("--haplotype", required=True, help="PanSN haplotype key, for the output")
    p.add_argument("--label", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-identity", type=float, default=0.90,
                   help="minimum gap-compressed identity for a hit to count (default 0.90)")
    p.add_argument("--min-frac", type=float, default=0.5,
                   help="merged aligned fraction at or above which a segment is called "
                        "NOT_PRIVATE (default 0.5)")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    stem = "%s.%s" % (a.label, a.haplotype.replace("#", "_").replace("/", "_"))

    qlen = {}
    # (query, target assembly) -> merged query intervals
    per_asm = defaultdict(list)
    best_id = defaultdict(float)
    n_rows = n_self = n_lowid = n_kept = 0
    seen_q = set()

    with popen(a.paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            n_rows += 1
            q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            tname = f[5]
            nmatch, blen = int(f[9]), int(f[10])
            qlen[q] = ql
            seen_q.add(q)

            asm = tname.split("::", 1)[0]
            if asm == a.self_assembly:
                n_self += 1
                continue
            ident = (nmatch / blen) if blen else 0.0
            if ident < a.min_identity:
                n_lowid += 1
                continue
            n_kept += 1
            per_asm[(q, asm)].append((qs, qe))
            if ident > best_id[q]:
                best_id[q] = ident

    # every segment in the FASTA should appear; those with no PAF row aligned nowhere
    fout = os.path.join(a.outdir, "%s.private_map.tsv" % stem)
    with open(fout, "w") as out:
        out.write("# per private segment: does this sequence exist in ANY OTHER assembly?\n")
        out.write("# Self-hits excluded by assembly tag (%s), not by identity -- a genuine\n"
                  % a.self_assembly)
        out.write("#   near-identical copy elsewhere is precisely the collapsed-repeat case\n")
        out.write("#   this is meant to detect, so an identity filter would hide it.\n")
        out.write("# aligned_frac is MERGED query coverage across all other assemblies;\n")
        out.write("#   summing PAF block lengths would double-count overlaps.\n")
        out.write("# verdict PRIVATE_CONFIRMED = no other assembly covers >= %.2f of it.\n"
                  % a.min_frac)
        out.write("haplotype\tsegment\tcontig\tstart\tend\tsegment_bp\t"
                  "n_other_assemblies\tbest_identity\taligned_frac_merged\t"
                  "max_frac_single_assembly\tverdict\n")

        rows = 0
        for q in sorted(seen_q):
            ql = qlen[q]
            # header is <haplotype>|<contig>|<start>|<end>|seg<N>
            parts = q.split("|")
            contig = parts[1] if len(parts) > 3 else "."
            start = parts[2] if len(parts) > 3 else "0"
            end = parts[3] if len(parts) > 3 else "0"
            seg = parts[4] if len(parts) > 4 else q

            asms = [k[1] for k in per_asm if k[0] == q]
            all_iv = []
            per_max = 0.0
            for asm in asms:
                iv = per_asm[(q, asm)]
                all_iv.extend(iv)
                fr = merged_len(iv) / ql if ql else 0.0
                per_max = max(per_max, fr)
            frac = (merged_len(all_iv) / ql) if ql else 0.0
            verdict = "NOT_PRIVATE" if frac >= a.min_frac else "PRIVATE_CONFIRMED"
            out.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%s\n"
                      % (a.haplotype, seg, contig, start, end, ql,
                         len(asms), best_id.get(q, 0.0), min(1.0, frac),
                         min(1.0, per_max), verdict))
            rows += 1

    with open(os.path.join(a.outdir, "%s.private_map_audit.tsv" % stem), "w") as out:
        out.write("metric\tvalue\n")
        out.write("haplotype\t%s\nself_assembly\t%s\n" % (a.haplotype, a.self_assembly))
        out.write("paf_rows\t%d\n" % n_rows)
        out.write("rows_self_excluded\t%d\n" % n_self)
        out.write("rows_below_identity_%.2f\t%d\n" % (a.min_identity, n_lowid))
        out.write("rows_kept\t%d\n" % n_kept)
        out.write("segments_with_any_paf_row\t%d\n" % len(seen_q))
        out.write("# segments_with_any_paf_row counts only segments minimap2 reported at all.\n")
        out.write("# Segments absent from the PAF aligned NOWHERE, including to their own\n")
        out.write("# assembly, which should not happen -- rows_self_excluded near zero is a\n")
        out.write("# warning that the --self-assembly tag does not match the index tags.\n")

    sys.stderr.write("[private_map] %s: %d PAF rows, %d self-excluded, %d kept, %d segments\n"
                     % (a.haplotype, n_rows, n_self, n_kept, len(seen_q)))
    if n_rows and n_self == 0:
        sys.stderr.write("WARNING: no self-hits excluded -- --self-assembly '%s' may not "
                         "match the index's assembly tags\n" % a.self_assembly)


if __name__ == "__main__":
    main()
