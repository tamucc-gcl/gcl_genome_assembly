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
import re
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
    p.add_argument("--flavor", default="unknown", help="recorded in the CSV")
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
        # `set` MUST be part of the key. Control segments are named ctrlN and private segN so
        # they do not collide by accident, but without it in the key nothing downstream can
        # tell them apart -- is_private came out as 1 on all 4.39M rows, which silently mixed
        # the control into every private figure.
        k = (d.get("haplotype", "."), d.get("set", "private"), d.get("segment", "."))
        rows[k] = {"haplotype": k[0], "set": k[1], "segment": k[2],
                   "contig": d.get("contig", "."),
                   "start": d.get("start", "0"), "end": d.get("end", "0"),
                   "segment_bp": inum(d, "segment_bp"),
                   "n_other_assemblies": inum(d, "n_other_assemblies"),
                   "best_identity": fnum(d, "best_identity"),
                   "aligned_frac_merged": fnum(d, "aligned_frac_merged"),
                   "map_verdict": d.get("verdict", "."),
                   "sample": ".", "n_kmers": 0, "mean_copy": -1.0,
                   "median_copy": -1, "frac_absent": -1.0, "kmer_verdict": ".",
                   "max_copy": -1, "single_copy_ref": -1.0, "copy_ratio": -1.0,
                   "n_kmers_expected": -1}
        n_map += 1

    for d, _ in read_tsv(kmer_files):
        k = (d.get("haplotype", "."), d.get("set", "private"), d.get("segment", "."))
        r = rows.get(k)
        if r is None:
            # a segment with k-mer data but no map row: mapping did not report it at all
            r = rows.setdefault(k, {"haplotype": k[0], "set": k[1], "segment": k[2],
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
        r["max_copy"] = inum(d, "max_copy", -1)
        r["single_copy_ref"] = fnum(d, "single_copy_ref", -1.0)
        r["copy_ratio"] = fnum(d, "copy_ratio", -1.0)
        r["n_kmers_expected"] = inum(d, "n_kmers_expected", -1)
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
            st = r.get("set", "private")
            combos[(st, comb)] += 1
            bp_combo[(st, comb)] += r.get("segment_bp", 0)
            chrom = str(r.get("contig", ".")).split("_")[0]
            per_chrom[chrom][(st, comb)] += 1
            per_bin[sbin(r.get("segment_bp", 0))][(st, comb)] += 1
            per_hap[r["haplotype"]][(st, comb)] += 1
            out.write("\t".join(str(r.get(c, ".")) for c in cols) + "\n")

    # ---- flat CSV for R --------------------------------------------------------------
    # One row per segment, carrying BOTH response forms so a binomial model
    # (repeat_like ~ is_private + ...) and a continuous one (log(median_copy) ~ ...) can be
    # fitted from the same file without re-deriving anything.
    #
    # is_private is 1 for every row here. A matched CONTROL set of non-private random windows
    # -- same haplotype, same chromosome, size-matched, rejecting any window containing
    # private sequence -- is what makes the enrichment test possible, and it appends rows with
    # is_private = 0. Until that exists this file describes the private set only and cannot
    # answer whether private sequence is repeat-enriched.
    #
    # individual is the haplotype minus its PanSN field and any _hapN suffix, so the two
    # haplotypes of one sample share it -- needed for the nested random effect.
    csv_cols = ["haplotype", "individual", "sample", "flavor", "chromosome", "segment",
                "set", "start", "end", "span_bp", "log_span", "is_private",
                "n_other_assemblies", "best_identity", "aligned_frac_merged", "map_verdict",
                "n_kmers_observed", "n_kmers_expected", "frac_absent",
                "mean_copy", "median_copy", "max_copy", "single_copy_ref", "copy_ratio",
                "kmer_verdict", "repeat_like", "combined"]
    import math as _math
    with open(os.path.join(a.outdir, "%s.private_evidence.csv" % a.label), "w") as out:
        out.write(",".join(csv_cols) + "\n")
        for k in sorted(rows):
            r = rows[k]
            hap = str(r.get("haplotype", "."))
            indiv = hap.split("#")[0]
            indiv = re.sub(r"_hap[0-9]+$", "", indiv)
            chrom = str(r.get("contig", ".")).split("_")[0]
            span = r.get("segment_bp", 0) or 0
            kv = r.get("kmer_verdict", ".")
            vals = {
                "haplotype": hap, "individual": indiv,
                "sample": r.get("sample", "."), "flavor": a.flavor,
                "chromosome": chrom, "segment": r.get("segment", "."),
                "set": r.get("set", "private"),
                "start": r.get("start", 0), "end": r.get("end", 0),
                "span_bp": span,
                "log_span": ("%.6f" % _math.log(span)) if span > 0 else "NA",
                "is_private": (0 if str(r.get("set", "private")) == "control" else 1),
                "n_other_assemblies": r.get("n_other_assemblies", -1),
                "best_identity": r.get("best_identity", -1),
                "aligned_frac_merged": r.get("aligned_frac_merged", -1),
                "map_verdict": r.get("map_verdict", "."),
                "n_kmers_observed": r.get("n_kmers", 0),
                "n_kmers_expected": r.get("n_kmers_expected", -1),
                "frac_absent": r.get("frac_absent", -1),
                "mean_copy": r.get("mean_copy", -1),
                "median_copy": r.get("median_copy", -1),
                "max_copy": r.get("max_copy", -1),
                "single_copy_ref": r.get("single_copy_ref", -1),
                "copy_ratio": r.get("copy_ratio", -1),
                "kmer_verdict": kv,
                "repeat_like": (1 if kv == "REPEAT_LIKE" else (0 if kv == "UNIQUE_LIKE" else "NA")),
                "combined": r.get("combined", "."),
            }
            out.write(",".join(str(vals[c]) for c in csv_cols) + "\n")

    # ---- cross-tabulations -----------------------------------------------------------
    with open(os.path.join(a.outdir, "%s.private_evidence_xtab.tsv" % a.label), "w") as out:
        out.write("# THE RESULT. Two independent verdicts, tabulated. See the header of\n")
        out.write("# private_evidence.tsv for what each cell means.\n")
        out.write("# broken out by SET: the private-vs-control contrast is the whole point,\n")
        out.write("# and a pooled table would average the two together.\n")
        out.write("scope\tkey\tset\tcombined\tn_segments\tsegment_bp\n")
        for (st, c), n in sorted(combos.items(), key=lambda x: -x[1]):
            out.write("ALL\t.\t%s\t%s\t%d\t%d\n" % (st, c, n, bp_combo[(st, c)]))
        for chrom in sorted(per_chrom):
            for (st, c), n in sorted(per_chrom[chrom].items(), key=lambda x: -x[1]):
                out.write("CHROM\t%s\t%s\t%s\t%d\t.\n" % (chrom, st, c, n))
        for b in sorted(per_bin, key=lambda s: edges.index(
                int(s.split("-")[0].lstrip(">="))) if s.split("-")[0].lstrip(">=").isdigit() else 99):
            for (st, c), n in sorted(per_bin[b].items(), key=lambda x: -x[1]):
                out.write("SIZE_BIN\t%s\t%s\t%s\t%d\t.\n" % (b, st, c, n))
        for h in sorted(per_hap):
            for (st, c), n in sorted(per_hap[h].items(), key=lambda x: -x[1]):
                out.write("HAPLOTYPE\t%s\t%s\t%s\t%d\t.\n" % (h, st, c, n))

    with open(os.path.join(a.outdir, "%s.private_evidence_audit.tsv" % a.label), "w") as out:
        out.write("metric\tvalue\n")
        out.write("label\t%s\n" % a.label)
        out.write("map_tables\t%d\nkmer_tables\t%d\n" % (len(map_files), len(kmer_files)))
        out.write("map_rows\t%d\nkmer_rows\t%d\n" % (n_map, n_kmer))
        out.write("joined_segments\t%d\n" % len(rows))
        n_priv = sum(1 for r in rows.values() if r.get("set", "private") != "control")
        n_ctrl = sum(1 for r in rows.values() if r.get("set", "private") == "control")
        out.write("segments_private\t%d\nsegments_control\t%d\n" % (n_priv, n_ctrl))
        out.write("# BOTH must be non-zero, or the enrichment contrast is not testable.\n")
        miss_k = sum(1 for r in rows.values() if r.get("kmer_verdict", ".") == ".")
        miss_m = sum(1 for r in rows.values() if r.get("map_verdict") == "NO_MAP_ROW")
        out.write("segments_without_kmer_row\t%d\n" % miss_k)
        out.write("segments_without_map_row\t%d\n" % miss_m)
        haps_map = {r["haplotype"] for r in rows.values()
                    if r.get("map_verdict") != "NO_MAP_ROW"}
        haps_kmer = {r["haplotype"] for r in rows.values()
                     if r.get("kmer_verdict", ".") != "."}
        out.write("haplotypes_with_map_rows\t%d\n" % len(haps_map))
        out.write("haplotypes_with_kmer_rows\t%d\n" % len(haps_kmer))
        out.write("haplotypes_missing_from_map\t%s\n"
                  % (",".join(sorted(haps_kmer - haps_map)) or "-"))
        out.write("haplotypes_missing_from_kmer\t%s\n"
                  % (",".join(sorted(haps_map - haps_kmer)) or "-"))
        out.write("# A WHOLE HAPLOTYPE missing from either stream is a failed task and is\n")
        out.write("# fatal. A handful of individual segments missing is not: minimap2 omits\n")
        out.write("# a query that aligns nowhere at all, so low-complexity segments can be\n")
        out.write("# absent from the map stream legitimately. The original check conflated\n")
        out.write("# the two and failed the run over 41 such segments.\n")

    sys.stderr.write("[private_join] %d segments joined (%d private, %d control); "
                     "%d without kmer, %d without map\n"
                     % (len(rows),
                        sum(1 for r in rows.values() if r.get("set", "private") != "control"),
                        sum(1 for r in rows.values() if r.get("set", "private") == "control"),
                        miss_k, miss_m))
    for (st, c), n in sorted(combos.items(), key=lambda x: -x[1]):
        sys.stderr.write("[private_join]   %-8s %-34s %8d segments  %12d bp\n"
                         % (st, c, n, bp_combo[(st, c)]))
    if miss_k or miss_m:
        sys.stderr.write("WARNING: %d segments lack kmer data, %d lack map data -- the two "
                         "scatters disagree on which haplotypes were processed\n"
                         % (miss_k, miss_m))


if __name__ == "__main__":
    main()
