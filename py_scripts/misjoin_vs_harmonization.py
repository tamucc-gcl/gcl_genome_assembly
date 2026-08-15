#!/usr/bin/env python3
"""
misjoin_vs_harmonization.py -- blind test of the Hi-C mis-join scan against an INDEPENDENT
predictor of where mis-joins are.

THE TEST
    The alignment-based harmonizer classifies every chromosome-scale scaffold without ever
    looking at Hi-C:

      1_single            one reference chromosome      -> no internal junction
      2_concordant        a fusion most voters agree on -> real junction, contiguous
      3_chimera_suspect   a fusion almost nobody else makes -> PREDICTED FALSE JOIN

    Group 3 is the positive set and group 1 is the negative set, chosen by a method that
    shares no data path with the Hi-C scan. If the scan finds deep localised dips in group 3
    and not in group 1, two independent methods agree on specific coordinates in assemblies
    neither was tuned on. Group 2 is the internal control: a real junction inside a
    composite should look like group 1, not group 3.

RAW RATIOS ARE NOT COMPARABLE ACROSS ASSEMBLIES
    Interior null median ratio, measured per assembly:

        CMat_203_hap1  0.922      CTlk_104_hap1  0.548
        CMat_203_hap2  0.918      CTlk_104_hap2  0.549
        CLim_110_hap1  0.942      CBau_104_hap1  0.638
        CLim_110_hap2  0.945

    That is the median over ~8000 interior points, so it is not mis-joins -- it is the
    expectation model. The representative separation for a distance bin is its midpoint
    (d = db*BIN + BIN/2), but the contact-distance decay is steep, so in the smallest bin
    the true mass sits well below the midpoint. That over-predicts the straddle window, and
    the steeper a library's decay the worse the over-prediction.

    Pooling raw ratios into one null would therefore make every CTlk scaffold look dipped
    against a null dominated by CMat and CLim points -- manufacturing exactly the false
    positives this test exists to find. So everything below uses the profile's
    ratio_over_scaffold_median column: a localised dip survives division by a median over
    ~800 points, while a whole-library offset cancels.

THE MINIMUM-STATISTIC PROBLEM
    A longer scaffold has more scan points, so its worst dip is deeper by chance. An
    earlier version corrected for this parametrically, p_min = 1 - (1 - F(x))^n, and that
    was worse than nothing: with few calibrating scaffolds F saturates, (1-F)^n underflows
    for large n, and on a synthetic fixture a genuinely CONTIGUOUS composite scored
    identically to a hard break -- the length correction reintroduced a length bias,
    inverted.

    Replaced with a like-for-like comparison: each scaffold's minimum is ranked against the
    distribution of MINIMA from 1_single scaffolds. Same selection effect on both sides, no
    parametric assumption. n_powered is reported alongside so residual length dependence
    stays visible rather than being silently "corrected".

    Also reported is a 3-point smoothed minimum. A genuine mis-join suppresses spanning
    contacts across the whole junction, so neighbours dip too; a lone deep point flanked by
    normal ones is noise.

Usage:
    misjoin_vs_harmonization.py --report 373251_harmonization_report.tsv \
        --profiles *.profile.tsv --out misjoin_validation.tsv
"""
import argparse
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, help="{taxid}.harmonization_report.tsv")
    ap.add_argument("--profiles", nargs="+", required=True)
    ap.add_argument("--min-expected", type=float, default=25.0)
    ap.add_argument("--min-scaffold-bp", type=int, default=8_000_000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # ---- harmonization classes (no Hi-C involved)
    grp, flags = {}, {}
    with open(a.report) as fh:
        hdr = None
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                I = {k: i for i, k in enumerate(hdr)}
                continue
            if int(f[I["length"]]) < a.min_scaffold_bp:
                continue
            fl = f[I["flags"]]
            cls = f[I["class"]]
            if cls == "chromosome":
                g = "1_single"
            elif "chimera_suspect" in fl:
                g = "3_chimera_suspect"
            elif "concordant" in fl:
                g = "2_concordant"
            else:
                g = "2_other_composite"
            key = (f[I["assembly"]], f[I["new_name"]])
            grp[key] = g
            flags[key] = fl

    # ---- profiles
    pts = defaultdict(list)
    for p in a.profiles:
        with open(p) as fh:
            hdr = None
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if hdr is None:
                    hdr = f
                    J = {k: i for i, k in enumerate(hdr)}
                    continue
                if float(f[J["expected"]]) < a.min_expected:
                    continue
                key = "ratio_over_scaffold_median"
                if key not in J:
                    sys.exit("profile %s lacks %s -- regenerate with the current "
                             "hic_misjoin_scan.py" % (p, key))
                v = f[J[key]]
                if v == "NA":
                    continue
                pts[(f[J["assembly"]], f[J["scaffold"]])].append(
                    (int(f[J["position"]]), float(v)))
    for k in pts:
        pts[k].sort()

    rows = []
    for k, v in pts.items():
        g = grp.get(k)
        if g is None or len(v) < 5:
            continue
        ratios = [r for _p, r in v]
        n = len(ratios)
        mi = min(range(n), key=lambda i: ratios[i])
        min_ratio, min_pos = ratios[mi], v[mi][0]
        sm = [sum(ratios[max(0, i - 1):i + 2]) / len(ratios[max(0, i - 1):i + 2])
              for i in range(n)]
        smi = min(range(n), key=lambda i: sm[i])
        rows.append({"asm": k[0], "scaf": k[1], "grp": g, "n": n,
                     "min_ratio": min_ratio, "min_pos": min_pos,
                     "smooth_min": sm[smi], "smooth_pos": v[smi][0],
                     "flags": flags.get(k, "-")})

    # like-for-like: rank each scaffold's minimum against the 1_single scaffolds' minima
    single_min = sorted(r["min_ratio"] for r in rows if r["grp"] == "1_single")
    single_smin = sorted(r["smooth_min"] for r in rows if r["grp"] == "1_single")
    if not single_min:
        sys.exit("no 1_single scaffolds in the profiles -- cannot calibrate")

    def pct_of(v, ref):
        return 100.0 * sum(1 for x in ref if x <= v) / float(len(ref))

    for r in rows:
        r["pct_min"] = pct_of(r["min_ratio"], single_min)
        r["pct_smooth"] = pct_of(r["smooth_min"], single_smin)

    rows.sort(key=lambda r: r["min_ratio"])
    with open(a.out, "w") as fh:
        fh.write("# values are ratio_over_scaffold_median, NOT raw ratio -- raw ratios "
                 "are not comparable across assemblies (null medians 0.55-0.95)\n")
        fh.write("# null of minima from %d 1_single scaffolds\n" % len(single_min))
        fh.write("assembly\tscaffold\tgroup\tn_powered\tmin_norm\tmin_position\t"
                 "pct_among_single_minima\tsmooth_min_norm\tsmooth_position\t"
                 "pct_smooth\tflags\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in (
                r["asm"], r["scaf"], r["grp"], r["n"], "%.4f" % r["min_ratio"],
                r["min_pos"], "%.1f" % r["pct_min"], "%.4f" % r["smooth_min"],
                r["smooth_pos"], "%.1f" % r["pct_smooth"], r["flags"])) + "\n")

    print("all values normalised per scaffold (ratio / that scaffold's median)")
    print("  %-20s %5s %10s %10s %10s %9s"
          % ("group", "n", "med min", "worst min", "med smooth", "med n_pts"))
    for g in sorted({r["grp"] for r in rows}):
        gr = [r for r in rows if r["grp"] == g]
        ms = sorted(r["min_ratio"] for r in gr)
        ss = sorted(r["smooth_min"] for r in gr)
        ns = sorted(r["n"] for r in gr)
        print("  %-20s %5d %10.3f %10.3f %10.3f %9d"
              % (g, len(gr), ms[len(ms) // 2], ms[0], ss[len(ss) // 2],
                 ns[len(ns) // 2]))
    print("\nDeepest 12 dips (pct = percentile among 1_single scaffolds' minima;")
    print("a mis-join should sit far below every single-chromosome scaffold):")
    print("  %-20s %-24s %-18s %8s %6s %8s %6s %7s"
          % ("assembly", "scaffold", "group", "min", "pct", "smooth", "pct", "n_pts"))
    for r in rows[:12]:
        print("  %-20s %-24s %-18s %8.4f %6.1f %8.4f %6.1f %7d"
              % (r["asm"], r["scaf"], r["grp"], r["min_ratio"], r["pct_min"],
                 r["smooth_min"], r["pct_smooth"], r["n"]))


if __name__ == "__main__":
    main()
