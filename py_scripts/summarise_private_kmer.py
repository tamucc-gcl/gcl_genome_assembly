#!/usr/bin/env python3
# ======================================================================================
# summarise_private_kmer.py
#
# Reads a `meryl-lookup -dump` stream and reports, per private segment, the copy number of
# its k-mers in the sample's own read database.
#
# WHAT THIS MEASURES AND WHY IT IS THE RIGHT MEASURE
# -------------------------------------------------
# The question is whether private sequence is repeat-derived. Copy number is taken from the
# sample's READ k-mer database, not from its assembly, deliberately: read multiplicity
# reflects true genomic copy number and is unaffected by whether the assembler collapsed the
# repeat. An assembly-derived count would be circular -- collapsed repeats are exactly what
# might be producing spurious private sequence in the first place.
#
# A segment whose k-mers sit at ~1x sample coverage is unique sequence. One sitting at many
# multiples is repeat-derived, and a private segment that is also high copy number is a
# strong candidate for graph collapse rather than biology -- which is the same conclusion
# PRIVATE_MAP reaches by a completely different route, so the two cross-check each other.
#
# This is a permanent part of the pipeline, not a placeholder. It needs no external
# annotation, no TE library and no reference database, so it works for any species the
# pipeline is pointed at. A RepeatMasker/RepeatModeler-based analysis can be layered on later
# from a separate annotation pipeline, but nothing here depends on that existing.
#
# WHY THE FORMAT IS AUTO-DETECTED
# -------------------------------
# meryl-lookup's -dump column layout is not something to hardcode from memory, and it has
# varied between meryl releases. The first data line is inspected: column 1 is taken as the
# sequence name and the last purely-numeric column as the k-mer value, with the detected
# layout written into the audit file so a wrong guess is visible rather than silent.
# --name-column / --value-column override it.
#
# Reads stdin so the dump is never written to disk -- roughly 137 Mb of sequence per
# haplotype means ~137M k-mer lines, which is minutes to stream but tens of GB to store.
#
# Dependencies: python3 stdlib.
#
# USAGE
#   meryl-lookup -dump -sequence x.fa -mers db.meryl \
#     | summarise_private_kmer.py --haplotype 'Sde-CBau_104#1' --sample Sde-CBau_104 \
#         --label 373251.full --outdir .
# ======================================================================================

import argparse
import os
import sys


def is_num(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--haplotype", required=True, help="PanSN haplotype key")
    p.add_argument("--sample", required=True, help="sample whose meryl DB was queried")
    p.add_argument("--label", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--name-column", type=int, default=0,
                   help="1-based; 0 = auto-detect (default)")
    p.add_argument("--value-column", type=int, default=0,
                   help="1-based; 0 = auto-detect as the last numeric column")
    p.add_argument("--repeat-threshold", type=float, default=3.0,
                   help="mean copy number at or above which a segment is flagged "
                        "REPEAT_LIKE (default 3.0)")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    stem = "%s.%s" % (a.label, a.haplotype.replace("#", "_").replace("/", "_"))

    # per segment: count, sum, and a value histogram for the median
    n_k, sum_k, zero_k = {}, {}, {}
    hist = {}
    name_col = a.name_column - 1 if a.name_column else None
    val_col = a.value_column - 1 if a.value_column else None
    detected = None
    n_lines = n_bad = 0

    for line in sys.stdin:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            f = line.split()
        if len(f) < 2:
            n_bad += 1
            continue

        if name_col is None or val_col is None:
            # auto-detect from the first usable line
            nc = 0
            vc = None
            for i in range(len(f) - 1, -1, -1):
                if is_num(f[i]):
                    vc = i
                    break
            if vc is None or vc == nc:
                n_bad += 1
                continue
            name_col, val_col = nc, vc
            detected = "%d columns; name=col%d value=col%d; first line: %s" % (
                len(f), nc + 1, vc + 1, "\t".join(f[:6]))
            sys.stderr.write("[private_kmer] detected %s\n" % detected)

        if val_col >= len(f) or name_col >= len(f):
            n_bad += 1
            continue
        try:
            v = float(f[val_col])
        except ValueError:
            n_bad += 1
            continue

        s = f[name_col]
        n_lines += 1
        n_k[s] = n_k.get(s, 0) + 1
        sum_k[s] = sum_k.get(s, 0.0) + v
        if v == 0:
            zero_k[s] = zero_k.get(s, 0) + 1
        # bounded histogram: exact for 0-99, then a tail bucket. Enough for a median and
        # avoids holding every value for 137M k-mers.
        h = hist.setdefault(s, [0] * 101)
        h[min(int(v), 100)] += 1

    def median_from_hist(h, total):
        half = total / 2.0
        run = 0
        for i, c in enumerate(h):
            run += c
            if run >= half:
                return i
        return 0

    fout = os.path.join(a.outdir, "%s.private_kmer.tsv" % stem)
    with open(fout, "w") as out:
        out.write("# per private segment: k-mer copy number in the SAMPLE'S OWN READ database.\n")
        out.write("# Read-derived, not assembly-derived, on purpose: read multiplicity reflects\n")
        out.write("#   true genomic copy number and is unaffected by assembly collapse, whereas\n")
        out.write("#   an assembly count would be circular -- collapse is one of the things\n")
        out.write("#   that may be producing spurious private sequence.\n")
        out.write("# median_copy is from a histogram capped at 100; values above that are\n")
        out.write("#   counted in the top bucket, so a median of 100 means '>=100'.\n")
        out.write("# frac_absent is k-mers with copy number 0 -- unexpected in a sample's own\n")
        out.write("#   reads, so a high value suggests the wrong meryl database was joined.\n")
        out.write("# REPEAT_LIKE at mean_copy >= %.1f. A segment that is BOTH private and\n"
                  % a.repeat_threshold)
        out.write("#   high copy number is a candidate for graph collapse rather than biology,\n")
        out.write("#   which PRIVATE_MAP tests independently -- the two should agree.\n")
        out.write("haplotype\tsample\tsegment\tcontig\tstart\tend\tn_kmers\t"
                  "mean_copy\tmedian_copy\tfrac_absent\tverdict\n")
        for s in sorted(n_k):
            n = n_k[s]
            mean = sum_k[s] / n if n else 0.0
            med = median_from_hist(hist[s], n)
            fa = zero_k.get(s, 0) / n if n else 0.0
            parts = s.split("|")
            contig = parts[1] if len(parts) > 3 else "."
            start = parts[2] if len(parts) > 3 else "0"
            end = parts[3] if len(parts) > 3 else "0"
            seg = parts[4] if len(parts) > 4 else s
            verdict = "REPEAT_LIKE" if mean >= a.repeat_threshold else "UNIQUE_LIKE"
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.4f\t%d\t%.4f\t%s\n"
                      % (a.haplotype, a.sample, seg, contig, start, end,
                         n, mean, med, fa, verdict))

    with open(os.path.join(a.outdir, "%s.private_kmer_audit.tsv" % stem), "w") as out:
        out.write("metric\tvalue\n")
        out.write("haplotype\t%s\nsample\t%s\n" % (a.haplotype, a.sample))
        out.write("kmer_lines\t%d\nunparseable_lines\t%d\n" % (n_lines, n_bad))
        out.write("segments\t%d\n" % len(n_k))
        out.write("detected_format\t%s\n" % (detected or "NONE -- no data parsed"))
        tot = sum(n_k.values()) or 1
        out.write("overall_frac_absent\t%.6f\n" % (sum(zero_k.values()) / tot))
        out.write("# detected_format records how the columns were interpreted. If it looks\n")
        out.write("# wrong, pass --name-column / --value-column rather than trusting this.\n")
        out.write("# overall_frac_absent should be near zero: these are the sample's own reads.\n")

    sys.stderr.write("[private_kmer] %s vs %s: %d k-mers over %d segments, %d unparseable\n"
                     % (a.haplotype, a.sample, n_lines, len(n_k), n_bad))
    if n_lines == 0:
        sys.stderr.write("ERROR: no k-mer lines parsed -- check meryl-lookup output format\n")
        sys.exit(1)
    if sum(zero_k.values()) / tot > 0.5:
        sys.stderr.write("WARNING: %.1f%% of k-mers absent from %s's own read database -- "
                         "the wrong meryl DB may have been joined\n"
                         % (100.0 * sum(zero_k.values()) / tot, a.sample))


if __name__ == "__main__":
    main()
