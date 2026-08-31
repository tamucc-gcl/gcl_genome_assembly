#!/usr/bin/env python3
# ======================================================================================
# join_private_evidence.py
#
# Merges the per-haplotype PRIVATE_MAP and PRIVATE_KMER tables for one graph flavour into a
# single per-segment table, and cross-tabulates their two verdicts.
#
# WHY ONE TABLE
# -------------
# The statistical tests on this are deliberately deferred, so the requirement is that when
# they happen they are a regression over columns that already exist rather than a fresh
# data-gathering exercise. Segment size, alignment evidence and copy number therefore land in
# one row, keyed on (haplotype, segment).
#
# WHY THE CROSS-TABULATION IS THE ACTUAL RESULT
# ---------------------------------------------
# PRIVATE_MAP and PRIVATE_KMER answer the same question by completely independent routes --
# one aligns the segment against the other assemblies, the other counts its k-mers in the
# sample's own reads. Neither sees the other's evidence, so their agreement is informative:
#
#   NOT_PRIVATE      x REPEAT_LIKE   the segment exists elsewhere AND is high copy number.
#                                    Graph collapse. Not biology.
#   PRIVATE_CONFIRMED x UNIQUE_LIKE  absent from other assemblies AND single copy in reads.
#                                    Genuinely novel sequence. The interesting class.
#   NOT_PRIVATE      x UNIQUE_LIKE   present elsewhere but single copy -- the graph failed to
#                                    merge homologous sequence. An alignment failure rather
#                                    than a repeat problem.
#   PRIVATE_CONFIRMED x REPEAT_LIKE  absent from other assemblies but high copy WITHIN this
#                                    sample. A haplotype-specific expansion, or a repeat family
#                                    the other assemblies collapsed away entirely.
#
# The off-diagonals are the ones worth looking at, and they are invisible unless the two are
# tabulated together.
#
# THE QUESTION THIS IS AIMED AT
# -----------------------------
# chr8 and chr9 carry 23.6% and 22.9% private sequence against an ~11% floor across the other
# chromosomes, consistently across all ten haplotypes and both graph flavours. Consistency
# across independently assembled haplotypes argues for biology; a per-chromosome breakdown of
# this cross-tabulation is what decides it.
#
# Dependencies: python3 stdlib.
#
# USAGE
#   join_private_evidence.py --map '*.private_map.tsv' --kmer '*.private_kmer.tsv' \
#       --label 373251.full --outdir .
# ======================================================================================

import argparse
import glob
import os
import sys
from collections import Counter, defaultdict


def read_tsv(paths):
    """Yield dict rows from '#'-commented TSVs, tolerating column drift."""
    for p in paths:
        hdr = None
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if hdr is None:
                    hdr = f
                    continue
                if len(f) != len(hdr):
                    continue
                yield dict(zip(hdr, f)), p


def expand(pats):
    out = []
    for p in pats:
        out.extend(sorted(glob.glob(p)) if any(c in p for c in "*?[") else [p])
    return [p for p in out if os.path.isfile(p)]


def fnum(d, k, dflt=0.0):
    try:
        return float(d.get(k, dflt))
    except (TypeError, ValueError):
        return dflt


def inum(d, k, dflt=0):
    try:
        return int(float(d.get(k, dflt)))
    except (TypeError, ValueError):
        return dflt


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--map", nargs="+", required=True)
    p.add_argument("--kmer", nargs="+", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--bins", default="1000,5000,10000,50000,100000,500000,1000000",
                   help="size bin edges, matching the SV and private-segment spectra so all "
                        "three are readable on one axis")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    map_files, kmer_files = expand(a.map), expand(a.kmer)
    sys.stderr.write("[private_join] %d map tables, %d kmer tables\n"
                     % (len(map_files), len(kmer_files)))

    edges = [0] + [int(x) for x in a.bins.split(",") if x.strip()] + [float("inf")]

    def sbin(x):
        for lo, hi in zip(edges, edges[1:]):
            if lo <= x < hi:
                return (">=%d" % lo) if hi == float("inf") else ("%d-%d" % (lo, hi))
        return "unbinned"

    rows = {}
    n_map = n_kmer = 0
    for d, _ in read_tsv(map_files):
        k = (d.get("haplotype", "."), d.get("segment", "."))
        rows[k] = {"haplotype": k[0], "segment": k[1],
                   "contig": d.get("contig", "."),
                   "start": d.get("start", "0"), "end": d.get("end", "0"),
                   "segment_bp": inum(d, "segment_bp"),
                   "n_other_assemblies": inum(d, "n_other_assemblies"),
                   "best_identity": fnum(d, "best_identity"),
                   "aligned_frac_merged": fnum(d, "aligned_frac_merged"),
                   "map_verdict": d.get("verdict", "."),
                   "sample": ".", "n_kmers": 0, "mean_copy": -1.0,
                   "median_copy": -1, "frac_absent": -1.0, "kmer_verdict": "."}
        n_map += 1

    for d, _ in read_tsv(kmer_files):
        k = (d.get("haplotype", "."), d.get("segment", "."))
        r = rows.get(k)
        if r is None:
            # a segment with k-mer data but no map row: mapping did not report it at all
            r = rows.setdefault(k, {"haplotype": k[0], "segment": k[1],
                                    "contig": d.get("contig", "."),
                                    "start": d.get("start", "0"), "end": d.get("end", "0"),
                                    "segment_bp": 0, "n_other_assemblies": -1,
                                    "best_identity": -1.0, "aligned_frac_merged": -1.0,
                                    "map_verdict": "NO_MAP_ROW"})
        r["sample"] = d.get("sample", ".")
        r["n_kmers"] = inum(d, "n_kmers")
        r["mean_copy"] = fnum(d, "mean_copy", -1.0)
        r["median_copy"] = inum(d, "median_copy", -1)
        r["frac_absent"] = fnum(d, "frac_absent", -1.0)
        r["kmer_verdict"] = d.get("verdict", ".")
        n_kmer += 1

    # ---- the joined per-segment table ------------------------------------------------
    cols = ["haplotype", "sample", "segment", "contig", "start", "end", "segment_bp",
            "n_other_assemblies", "best_identity", "aligned_frac_merged", "map_verdict",
            "n_kmers", "mean_copy", "median_copy", "frac_absent", "kmer_verdict",
            "combined"]
    combos = Counter()
    per_chrom = defaultdict(Counter)
    per_bin = defaultdict(Counter)
    per_hap = defaultdict(Counter)
    bp_combo = Counter()

    with open(os.path.join(a.outdir, "%s.private_evidence.tsv" % a.label), "w") as out:
        out.write("# one row per private segment, both evidence streams joined.\n")
        out.write("# map_verdict comes from aligning the segment against the OTHER assemblies;\n")
        out.write("# kmer_verdict from its copy number in the sample's OWN reads. Independent.\n")
        out.write("# combined is the pair, and its off-diagonals are the informative cases:\n")
        out.write("#   NOT_PRIVATE+REPEAT_LIKE        graph collapse, not biology\n")
        out.write("#   PRIVATE_CONFIRMED+UNIQUE_LIKE  genuinely novel sequence\n")
        out.write("#   NOT_PRIVATE+UNIQUE_LIKE        alignment failure, not a repeat problem\n")
        out.write("#   PRIVATE_CONFIRMED+REPEAT_LIKE  haplotype-specific expansion\n")
        out.write("# -1 means the stream had no row for that segment.\n")
        out.write("\t".join(cols) + "\n")
        for k in sorted(rows):
            r = rows[k]
            mv, kv = r.get("map_verdict", "."), r.get("kmer_verdict", ".")
            comb = "%s+%s" % (mv, kv)
            r["combined"] = comb
            combos[comb] += 1
            bp_combo[comb] += r.get("segment_bp", 0)
            chrom = str(r.get("contig", ".")).split("_")[0]
            per_chrom[chrom][comb] += 1
            per_bin[sbin(r.get("segment_bp", 0))][comb] += 1
            per_hap[r["haplotype"]][comb] += 1
            out.write("\t".join(str(r.get(c, ".")) for c in cols) + "\n")

    # ---- cross-tabulations -----------------------------------------------------------
    with open(os.path.join(a.outdir, "%s.private_evidence_xtab.tsv" % a.label), "w") as out:
        out.write("# THE RESULT. Two independent verdicts, tabulated. See the header of\n")
        out.write("# private_evidence.tsv for what each cell means.\n")
        out.write("scope\tkey\tcombined\tn_segments\tsegment_bp\n")
        for c, n in sorted(combos.items(), key=lambda x: -x[1]):
            out.write("ALL\t.\t%s\t%d\t%d\n" % (c, n, bp_combo[c]))
        for chrom in sorted(per_chrom):
            for c, n in sorted(per_chrom[chrom].items(), key=lambda x: -x[1]):
                out.write("CHROM\t%s\t%s\t%d\t.\n" % (chrom, c, n))
        for b in sorted(per_bin, key=lambda s: edges.index(
                int(s.split("-")[0].lstrip(">="))) if s.split("-")[0].lstrip(">=").isdigit() else 99):
            for c, n in sorted(per_bin[b].items(), key=lambda x: -x[1]):
                out.write("SIZE_BIN\t%s\t%s\t%d\t.\n" % (b, c, n))
        for h in sorted(per_hap):
            for c, n in sorted(per_hap[h].items(), key=lambda x: -x[1]):
                out.write("HAPLOTYPE\t%s\t%s\t%d\t.\n" % (h, c, n))

    with open(os.path.join(a.outdir, "%s.private_evidence_audit.tsv" % a.label), "w") as out:
        out.write("metric\tvalue\n")
        out.write("label\t%s\n" % a.label)
        out.write("map_tables\t%d\nkmer_tables\t%d\n" % (len(map_files), len(kmer_files)))
        out.write("map_rows\t%d\nkmer_rows\t%d\n" % (n_map, n_kmer))
        out.write("joined_segments\t%d\n" % len(rows))
        miss_k = sum(1 for r in rows.values() if r.get("kmer_verdict", ".") == ".")
        miss_m = sum(1 for r in rows.values() if r.get("map_verdict") == "NO_MAP_ROW")
        out.write("segments_without_kmer_row\t%d\n" % miss_k)
        out.write("segments_without_map_row\t%d\n" % miss_m)
        out.write("# Both should be 0. The two streams run on the SAME FASTAs, so a segment\n")
        out.write("# present in one and absent from the other means a task failed or a\n")
        out.write("# haplotype was dropped from one scatter but not the other.\n")

    sys.stderr.write("[private_join] %d segments joined; %d without kmer, %d without map\n"
                     % (len(rows), miss_k, miss_m))
    for c, n in sorted(combos.items(), key=lambda x: -x[1]):
        sys.stderr.write("[private_join]   %-40s %8d segments  %12d bp\n" % (c, n, bp_combo[c]))
    if miss_k or miss_m:
        sys.stderr.write("WARNING: %d segments lack kmer data, %d lack map data -- the two "
                         "scatters disagree on which haplotypes were processed\n"
                         % (miss_k, miss_m))


if __name__ == "__main__":
    main()
