#!/usr/bin/env python3
# ======================================================================================
# summarise_private_kmer.py
#
# Reads `meryl-lookup -wig-count` and reports, per private segment, k-mer copy number in the
# sample's own read database -- expressed as a MULTIPLE of that haplotype's single-copy
# coverage, not as a raw count.
#
# WHY -wig-count
# --------------
# meryl-lookup's report types are -bed, -bed-runs, -wig-count, -wig-depth, -existence,
# -include and -exclude. There is no -dump; an earlier version of this script invented one
# from memory and every task failed. -wig-count is the right one: "the multiplicity of the
# kmer starting at each position in the sequence, if it exists in an input kmer database".
# -existence gives only present/absent counts, which is not a copy-number measure.
#
# FORMAT (from `meryl-lookup -wig-count -help`, verified against real output)
#     variableStep chrom=<sequence_name>
#     <position>\t<multiplicity>
# Writes to stdout when -output is omitted, so this can be streamed. The chrom= line carries
# the FULL FASTA header, which extract_private_fasta.py builds as
#     <haplotype>|<contig>|<start>|<end>|seg<N>
# so segment identity, coordinates and length are all available without a lookup.
#
# POSITIONS WITH NO DATABASE HIT ARE OMITTED, not reported as zero. So absence is measured as
# (expected k-mer positions) - (lines seen), where expected = (end - start) - k + 1 from the
# header coordinates. Counting only the lines present would report frac_absent = 0 always.
#
# THE THRESHOLD IS SELF-CALIBRATING, AND HAS TO BE
# ------------------------------------------------
# Copy number here is READ multiplicity, so its scale is sequencing depth. Measured on
# Sde-CBau_104: single-copy k-mers sit at ~14, and a satellite segment at chr10 position 0
# reads ~360,000. An absolute threshold is therefore meaningless -- the original 3.0 would
# have called every segment REPEAT_LIKE at 14x coverage -- and it cannot be a shared constant
# either, because depth differs per sample.
#
# So the reference level is taken from the DATA: the median of per-segment median
# multiplicities across this haplotype, which is the single-copy mode for anything but a
# genome that is mostly repeat. copy_ratio = segment median / that reference, and
# REPEAT_LIKE is copy_ratio >= --repeat-ratio (default 3.0, now a MULTIPLE rather than a
# count). The reference is reported so the calibration is auditable.
#
# Two passes are unavoidable for that, but the second is over ~170k per-segment summaries
# rather than the ~137M k-mer lines, so the stream is still read once.
#
# Dependencies: python3 stdlib.
#
# USAGE
#   meryl-lookup -wig-count -sequence x.private.fa -mers db.meryl \
#     | summarise_private_kmer.py --haplotype 'Sde-CBau_104#1' --sample Sde-CBau_104 \
#         --label 373251.full --outdir . --kmer 21
# ======================================================================================

import argparse
import os
import sys


def median(sorted_vals):
    n = len(sorted_vals)
    if n == 0:
        return 0.0
    mid = n // 2
    if n % 2:
        return float(sorted_vals[mid])
    return (sorted_vals[mid - 1] + sorted_vals[mid]) / 2.0


def parse_header(chrom):
    """<haplotype>|<contig>|<start>|<end>|seg<N> -> (seg, contig, start, end, span)."""
    p = chrom.split("|")
    if len(p) >= 5:
        try:
            start, end = int(p[2]), int(p[3])
        except ValueError:
            start, end = 0, 0
        return p[4], p[1], start, end, max(0, end - start)
    return chrom, ".", 0, 0, 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--haplotype", required=True, help="PanSN haplotype key")
    ap.add_argument("--sample", required=True, help="sample whose meryl DB was queried")
    ap.add_argument("--label", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--kmer", type=int, default=21,
                    help="k used to build the meryl database (default 21). Only affects the "
                         "expected k-mer count per segment, hence frac_absent.")
    ap.add_argument("--repeat-ratio", type=float, default=3.0,
                    help="copy number relative to this haplotype's single-copy level at or "
                         "above which a segment is REPEAT_LIKE (default 3.0). A MULTIPLE, "
                         "not a raw count: measured single-copy depth was ~14x on this "
                         "cohort, so an absolute threshold of 3 would flag everything.")
    ap.add_argument("--min-kmers", type=int, default=10,
                    help="segments with fewer observed k-mers than this are reported but "
                         "excluded from the single-copy reference (default 10)")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    stem = "%s.%s" % (a.label, a.haplotype.replace("#", "_").replace("/", "_"))

    # ---- pass 1: stream the wiggle, one summary per segment --------------------------
    segs = []            # (seg, contig, start, end, span, n_obs, total, med, mx)
    cur = None
    vals = []
    n_lines = n_bad = n_headers = 0

    def close_segment():
        if cur is None:
            return
        seg, contig, start, end, span = cur
        vals.sort()
        segs.append((seg, contig, start, end, span, len(vals),
                     sum(vals), median(vals), (vals[-1] if vals else 0)))

    for line in sys.stdin:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith("variableStep"):
            close_segment()
            n_headers += 1
            chrom = ""
            for tok in line.split():
                if tok.startswith("chrom="):
                    chrom = tok[6:]
            # a chrom= value can itself contain '=' or spaces in principle; take everything
            # after the first 'chrom=' to be safe
            i = line.find("chrom=")
            if i >= 0:
                chrom = line[i + 6:].strip()
            cur = parse_header(chrom)
            vals = []
            continue
        f = line.split("\t")
        if len(f) < 2:
            f = line.split()
        if len(f) < 2:
            n_bad += 1
            continue
        try:
            vals.append(int(f[1]))
        except ValueError:
            n_bad += 1
            continue
        n_lines += 1
    close_segment()

    if n_headers == 0:
        sys.stderr.write("ERROR: no 'variableStep' header seen. meryl-lookup -wig-count "
                         "should emit one per sequence; check the report type and that the "
                         "FASTA was not empty.\n")
        sys.exit(1)

    # ---- single-copy reference, from the data ----------------------------------------
    # K-MER-WEIGHTED median of per-segment medians, not the plain median: a segment
    # contributes in proportion to how much sequence it represents. The unweighted version
    # lets a 50-kmer fragment count as much as a 3,000-kmer one, which on a small segment set
    # lands the reference between modes -- a 4-segment fixture with 14x and 44x segments gave
    # 29.5, calibrating everything wrongly.
    #
    # Still a median rather than a mean, so a minority of extreme-copy segments cannot move
    # it: one 390 kb satellite at ~360,000 must not become the single-copy level.
    wq = sorted(((s[7], s[5]) for s in segs if s[5] >= a.min_kmers), key=lambda x: x[0])
    tot_w = sum(w for _, w in wq)
    ref_level = 0.0
    if tot_w:
        half, run = tot_w / 2.0, 0
        for val, w in wq:
            run += w
            if run >= half:
                ref_level = float(val)
                break
    per_seg_med = [v for v, _ in wq]
    if ref_level <= 0:
        sys.stderr.write("WARNING: single-copy reference computed as %.3f; copy_ratio will be "
                         "reported as -1 and no REPEAT_LIKE calls made\n" % ref_level)

    # ---- pass 2: write, normalising against the reference ----------------------------
    fout = os.path.join(a.outdir, "%s.private_kmer.tsv" % stem)
    n_repeat = n_unique = 0
    with open(fout, "w") as out:
        out.write("# per private segment: k-mer copy number in the SAMPLE'S OWN READ database.\n")
        out.write("# Read-derived, not assembly-derived, deliberately: read multiplicity\n")
        out.write("#   reflects true genomic copy number and is unaffected by assembly\n")
        out.write("#   collapse, whereas an assembly count would be circular -- collapse is\n")
        out.write("#   one of the things that may be producing spurious private sequence.\n")
        out.write("#\n")
        out.write("# copy_ratio is median_copy / single_copy_reference, where the reference is\n")
        out.write("#   the median of per-segment medians for THIS haplotype (%.3f here).\n"
                  % ref_level)
        out.write("#   Raw multiplicity is sequencing depth and so is not comparable between\n")
        out.write("#   samples: single-copy measured ~14x on this cohort while a satellite\n")
        out.write("#   segment read ~360,000. REPEAT_LIKE is copy_ratio >= %.2f.\n"
                  % a.repeat_ratio)
        out.write("#\n")
        out.write("# frac_absent uses EXPECTED k-mer positions (span - k + 1) because\n")
        out.write("#   -wig-count OMITS positions with no database hit rather than\n")
        out.write("#   reporting zero. Counting only observed lines would always give 0.\n")
        out.write("#   High absence means the wrong meryl database was joined: these are the\n")
        out.write("#   sample's own reads, so its own k-mers must be present.\n")
        out.write("haplotype\tsample\tsegment\tcontig\tstart\tend\tspan_bp\t"
                  "n_kmers_observed\tn_kmers_expected\tfrac_absent\t"
                  "mean_copy\tmedian_copy\tmax_copy\tsingle_copy_ref\tcopy_ratio\tverdict\n")
        tot_obs = tot_exp = 0
        for seg, contig, start, end, span, n_obs, total, med, mx in segs:
            expected = max(0, span - a.kmer + 1)
            tot_obs += n_obs
            tot_exp += expected
            frac_abs = (1.0 - (n_obs / expected)) if expected else -1.0
            frac_abs = max(0.0, min(1.0, frac_abs)) if expected else -1.0
            mean = (total / n_obs) if n_obs else 0.0
            ratio = (med / ref_level) if ref_level > 0 else -1.0
            if ratio < 0:
                verdict = "UNCALIBRATED"
            elif ratio >= a.repeat_ratio:
                verdict = "REPEAT_LIKE"
                n_repeat += 1
            else:
                verdict = "UNIQUE_LIKE"
                n_unique += 1
            out.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%.4f\t%.1f\t%d\t%.4f\t%s\t%s\n"
                      % (a.haplotype, a.sample, seg, contig, start, end, span,
                         n_obs, expected,
                         ("%.4f" % frac_abs) if frac_abs >= 0 else "NA",
                         mean, med, mx, ref_level,
                         ("%.4f" % ratio) if ratio >= 0 else "NA", verdict))

    overall_abs = (1.0 - (tot_obs / tot_exp)) if tot_exp else -1.0
    with open(os.path.join(a.outdir, "%s.private_kmer_audit.tsv" % stem), "w") as out:
        out.write("metric\tvalue\n")
        out.write("haplotype\t%s\nsample\t%s\n" % (a.haplotype, a.sample))
        out.write("kmer\t%d\n" % a.kmer)
        out.write("wiggle_headers\t%d\n" % n_headers)
        out.write("kmer_lines\t%d\nunparseable_lines\t%d\n" % (n_lines, n_bad))
        out.write("segments\t%d\n" % len(segs))
        out.write("single_copy_reference\t%.4f\n" % ref_level)
        out.write("segments_in_reference\t%d\n" % len(per_seg_med))
        out.write("kmers_in_reference\t%d\n" % tot_w)
        out.write("# the reference is the KMER-WEIGHTED median of per-segment medians, so a\n")
        out.write("# segment counts in proportion to the sequence it represents.\n")
        out.write("repeat_ratio_threshold\t%.2f\n" % a.repeat_ratio)
        out.write("verdict_REPEAT_LIKE\t%d\nverdict_UNIQUE_LIKE\t%d\n" % (n_repeat, n_unique))
        out.write("kmers_observed\t%d\nkmers_expected\t%d\n" % (tot_obs, tot_exp))
        out.write("overall_frac_absent\t%s\n"
                  % (("%.6f" % max(0.0, overall_abs)) if overall_abs >= 0 else "NA"))
        out.write("# single_copy_reference is the median of per-segment medians and is the\n")
        out.write("# calibration for copy_ratio. If it is far from the sample's expected read\n")
        out.write("# depth, the meryl database or the join is wrong.\n")
        out.write("# overall_frac_absent should be near zero -- the sample's own reads.\n")

    sys.stderr.write("[private_kmer] %s vs %s: %d segments, %d k-mer lines, %d unparseable\n"
                     % (a.haplotype, a.sample, len(segs), n_lines, n_bad))
    sys.stderr.write("[private_kmer]   single-copy reference %.2f (kmer-weighted over %d "
                     "segments / %d kmers)\n" % (ref_level, len(per_seg_med), tot_w))
    sys.stderr.write("[private_kmer]   REPEAT_LIKE %d / UNIQUE_LIKE %d at ratio >= %.2f\n"
                     % (n_repeat, n_unique, a.repeat_ratio))
    if overall_abs > 0.5:
        sys.stderr.write("WARNING: %.1f%% of expected k-mers absent from %s's own read "
                         "database -- check the sample join\n"
                         % (100.0 * overall_abs, a.sample))


if __name__ == "__main__":
    main()
