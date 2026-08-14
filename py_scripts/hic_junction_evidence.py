#!/usr/bin/env python3
"""
hic_junction_evidence.py -- does Hi-C say two scaffolds of ONE assembly are a single
chromosome the scaffolder failed to join?

A WITHIN-assembly measurement, independent of the harmonization frame, of the reference
choice, and of every other assembly -- which is what lets it speak to junctions the
cross-assembly alignment consensus cannot decide.

--------------------------------------------------------------------------------------
THE NULL: OBSERVED COVERAGE, NOT AREA
--------------------------------------------------------------------------------------
Expected contacts in a corner are proportional to the coverage actually present in those
windows, not to their area:

    expected = total_AB * (cov_A[window] / cov_A[total]) * (cov_B[window] / cov_B[total])

This is not a refinement -- an area-based null is unusable on real scaffold ends. Measured
on Sde-CBau_104_hap2, contact density relative to each chromosome's own mean runs:

    0-100 kb from the end   0.31      <- ~70% depleted
    100-200 kb              0.71
    200-300 kb              0.92
    300-400 kb              1.05      <- recovered
    1-5 Mb                  1.4-1.9   <- sub-terminal EXCESS

A ~6x swing along a chromosome. Under an area null a 200 kb terminal window predicted 6.6
contacts where coverage supports ~0.6, so observing zero looked like extreme depletion
(enrichment 0.0000) when it was unremarkable. At the other extreme a 10 Mb window sits
mostly inside the 1.4-1.9x zone, so everything looked mildly enriched -- which is why the
first run's apparent 1.22-1.25 "signal" on the three known junctions cannot be trusted: it
is the same magnitude as the coverage artifact it was measured against.

Conditioning on marginal coverage removes both failure modes. A window with no coverage
predicts no contacts and contributes no evidence, rather than manufacturing depletion.

--------------------------------------------------------------------------------------
THE POSITIVE CONTROL: ARTIFICIAL SPLITS
--------------------------------------------------------------------------------------
A null says what a non-junction looks like. It does not say whether this library, at this
depth, with this window, can detect a junction at all. Without that, "no support" is
uninterpretable -- and it is a live outcome here, because YaHS had this same Hi-C and
declined to make these joins.

So each large scaffold is split at its midpoint and the two halves scored with the identical
statistic. Those halves are a KNOWN single chromosome, so they give the empirical
distribution of true-junction signal in this dataset. Candidates are then read against two
distributions:

    inter-scaffold pairs   -> the null      (what "separate" looks like)
    artificial splits      -> the positive  (what "joined" looks like)

Interpretation follows without tuning:

    splits high, candidate high  -> real junction
    splits high, candidate low   -> Hi-C refutes the join
    splits low                   -> the test is insensitive here. Say so and stop. Do NOT
                                    lower thresholds until something passes.

The control is free: split halves are scored from the same streaming pass, using the
intra-scaffold contacts the junction test discards anyway.

--------------------------------------------------------------------------------------
WINDOWS
--------------------------------------------------------------------------------------
The window is [--window-inner, --window-inner + --window-bp) from each scaffold end. The
inner offset skips the dead zone. Both are recorded in the output; --scan sweeps them and
makes no calls.

Usage:
    hic_junction_evidence.py --pairs in.pairs.gz --fai asm.fa.fai \
        --assembly-id ID --out ID.junction_evidence.tsv
    hic_junction_evidence.py ... --scan          # calibration sweep, no calls
"""
import argparse
import gzip
import sys
from collections import defaultdict

CORNERS = (("start", "start"), ("start", "end"), ("end", "start"), ("end", "end"))
COVBIN = 100_000


def openf(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def quantiles(vals, qs):
    if not vals:
        return [0.0] * len(qs)
    s = sorted(vals)
    out = []
    for q in qs:
        k = (len(s) - 1) * q
        lo = int(k)
        hi = min(lo + 1, len(s) - 1)
        out.append(s[lo] + (s[hi] - s[lo]) * (k - lo))
    return out


def auc(pos_vals, neg_vals):
    """P(control > null). Scale-free, cannot blow up, reported as a diagnostic only.

    NOT used to set a threshold. Two reasons, both learned the hard way here:

      * The null CONTAINS the junctions being searched for. Any threshold chosen to give
        zero false positives therefore excludes the target by construction -- a Youden
        threshold on this data sat above the planted junction in testing.
      * An artificial mid-chromosome split is a STRONGER positive than a real missed
        junction: the two halves of an intact chromosome are in maximal proximity, whereas
        two scaffolds that failed to join may be separated by unassembled sequence or a
        repeat that suppresses contact. So the controls bound what is achievable, they do
        not estimate what a real junction scores.

    Thresholds come from --min-z (a calibrated Poisson deviate, data-independent) and
    --min-enrichment. The controls measure SENSITIVITY against those fixed thresholds; the
    null shows what clears them. Nothing is tuned to the data.
    """
    if not pos_vals or not neg_vals:
        return 0.0
    wins = ties = 0
    for p in pos_vals:
        for n in neg_vals:
            if p > n:
                wins += 1
            elif p == n:
                ties += 1
    return (wins + 0.5 * ties) / float(len(pos_vals) * len(neg_vals))


def window_cov(covbins, lo, hi):
    """Coverage (contact ends) in [lo, hi) from COVBIN-sized marginal bins."""
    if hi <= lo:
        return 0
    tot = 0
    for b in range(lo // COVBIN, (hi - 1) // COVBIN + 1):
        tot += covbins.get(b, 0)
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--assembly-id", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-scaffold-bp", type=int, default=1_000_000)
    ap.add_argument("--window-bp", type=int, default=1_000_000)
    ap.add_argument("--window-inner", type=int, default=300_000,
                    help="skip this much from the scaffold end before the window starts, "
                         "covering the terminal coverage dead zone")
    ap.add_argument("--min-total-contacts", type=int, default=100)
    ap.add_argument("--min-window-cov", type=int, default=200,
                    help="a corner whose windows carry less coverage than this is reported "
                         "low_coverage and never called; absence of contact is not "
                         "evidence of separation")
    ap.add_argument("--max-candidates", type=int, default=1000)
    ap.add_argument("--min-enrichment", type=float, default=2.0)
    ap.add_argument("--min-z", type=float, default=5.0)
    ap.add_argument("--control-min-scaffold-bp", type=int, default=8_000_000)
    ap.add_argument("--scan", action="store_true")
    a = ap.parse_args()

    lens = {}
    with open(a.fai) as fh:
        for line in fh:
            if line.strip():
                f = line.split("\t")
                lens[f[0]] = int(f[1])
    cand = sorted((n for n, L in lens.items() if L >= a.min_scaffold_bp),
                  key=lambda n: (-lens[n], n))
    if len(cand) > a.max_candidates:
        sys.stderr.write("[hic-junction %s] %d candidates trimmed to %d\n"
                         % (a.assembly_id, len(cand), a.max_candidates))
        cand = cand[:a.max_candidates]
    idx = {n: i for i, n in enumerate(cand)}
    ctrl = [n for n in cand if lens[n] >= a.control_min_scaffold_bp]
    mid = {n: lens[n] // 2 for n in ctrl}

    cov = defaultdict(lambda: defaultdict(int))
    covtot = defaultdict(int)
    total = defaultdict(int)
    pos = defaultdict(list)
    ctrl_total = defaultdict(int)
    ctrl_pos = defaultdict(list)
    keep = 3_000_000 if a.scan else a.window_inner + a.window_bp
    n_pairs = 0

    with openf(a.pairs) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_pairs += 1
            f = line.rstrip("\n").split("\t")
            c1, c2 = f[1], f[3]
            i1, i2 = idx.get(c1), idx.get(c2)
            if i1 is None or i2 is None:
                continue
            p1, p2 = int(f[2]), int(f[4])
            cov[c1][(p1 - 1) // COVBIN] += 1
            cov[c2][(p2 - 1) // COVBIN] += 1
            covtot[c1] += 1
            covtot[c2] += 1
            if c1 == c2:
                m = mid.get(c1)
                if m is not None and ((p1 <= m) != (p2 <= m)):
                    ctrl_total[c1] += 1
                    d1, d2 = abs(p1 - m), abs(p2 - m)
                    if d1 < keep and d2 < keep:
                        ctrl_pos[c1].append((d1, d2))
                continue
            if i1 <= i2:
                ka, kb, pa, pb = c1, c2, p1, p2
            else:
                ka, kb, pa, pb = c2, c1, p2, p1
            total[(ka, kb)] += 1
            da = min(pa - 1, lens[ka] - pa)
            db = min(pb - 1, lens[kb] - pb)
            if da < keep and db < keep:
                pos[(ka, kb)].append((pa, pb))

    def score_pair(ka, kb, inner, width):
        La, Lb = lens[ka], lens[kb]
        tot = total[(ka, kb)]
        if tot == 0:
            return None
        wa = {"start": (inner, inner + width), "end": (La - inner - width, La - inner)}
        wb = {"start": (inner, inner + width), "end": (Lb - inner - width, Lb - inner)}
        fa = {e: window_cov(cov[ka], max(0, lo), min(La, hi)) / max(1, covtot[ka])
              for e, (lo, hi) in wa.items()}
        fb = {e: window_cov(cov[kb], max(0, lo), min(Lb, hi)) / max(1, covtot[kb])
              for e, (lo, hi) in wb.items()}
        cnt = defaultdict(int)
        for pa, pb in pos[(ka, kb)]:
            ea = ("start" if wa["start"][0] <= pa - 1 < wa["start"][1] else
                  ("end" if wa["end"][0] <= pa - 1 < wa["end"][1] else None))
            if ea is None:
                continue
            eb = ("start" if wb["start"][0] <= pb - 1 < wb["start"][1] else
                  ("end" if wb["end"][0] <= pb - 1 < wb["end"][1] else None))
            if eb is None:
                continue
            cnt[(ea, eb)] += 1
        best = None
        for c in CORNERS:
            exp = tot * fa[c[0]] * fb[c[1]]
            if exp <= 0:
                continue
            obs = cnt.get(c, 0)
            enr = obs / exp
            if best is None or enr > best["enr"]:
                best = {"corner": c, "obs": obs, "exp": exp, "enr": enr,
                        "z": (obs - exp) / (exp ** 0.5), "tot": tot,
                        "wcov": min(fa[c[0]] * covtot[ka], fb[c[1]] * covtot[kb])}
        return best

    def score_control(c, inner, width):
        L = lens[c]
        tot = ctrl_total[c]
        if tot == 0:
            return None
        m = mid[c]
        lo, hi = inner, inner + width
        cA = window_cov(cov[c], max(0, m - hi), max(0, m - lo))
        cB = window_cov(cov[c], min(L, m + lo), min(L, m + hi))
        totA = window_cov(cov[c], 0, m)
        totB = window_cov(cov[c], m, L)
        if totA <= 0 or totB <= 0:
            return None
        exp = tot * (cA / totA) * (cB / totB)
        if exp <= 0:
            return None
        obs = sum(1 for d1, d2 in ctrl_pos[c] if lo <= d1 < hi and lo <= d2 < hi)
        return {"obs": obs, "exp": exp, "enr": obs / exp,
                "z": (obs - exp) / (exp ** 0.5), "tot": tot, "wcov": min(cA, cB)}

    if a.scan:
        def collect(inner, width):
            cs, ns = [], []
            for c in ctrl:
                r = score_control(c, inner, width)
                if r and r["wcov"] >= a.min_window_cov:
                    cs.append((r["enr"], r["z"], c))
            for i, ka in enumerate(cand):
                for kb in cand[i + 1:]:
                    r = score_pair(ka, kb, inner, width)
                    if (r and r["tot"] >= a.min_total_contacts
                            and r["wcov"] >= a.min_window_cov):
                        ns.append((r["enr"], r["z"], "%s+%s" % (ka, kb)))
            return cs, ns

        print("CONTROLS are artificial mid-chromosome splits -- KNOWN junctions in this")
        print("library. 'thr' is the control p10, i.e. the threshold that keeps 90%% of")
        print("them. n_hi counts inter-scaffold pairs reaching it; n_hiz adds z >= %.1f."
              % a.min_z)
        print("Those are CANDIDATES, not false positives -- named under the table.")
        print("A small n_hiz with a high thr is what a usable setting looks like.\n")
        hdr = ("%-7s %-8s | %7s %7s %7s | %7s %7s | %5s %8s %6s"
               % ("inner", "width", "c_med", "c_medz", "sens", "n_med", "n_medz",
                  "AUC", "n_pass", "n_ctrl"))
        print(hdr)
        print("-" * len(hdr))
        best = None
        for inner in (0, 100_000, 300_000):
            for width in (50_000, 100_000, 250_000, 500_000, 1_000_000, 2_000_000):
                cs, ns = collect(inner, width)
                if not cs or not ns:
                    continue
                def passes(e, z):
                    return e >= a.min_enrichment and z >= a.min_z
                c_med = quantiles([x[0] for x in cs], [0.5])[0]
                c_medz = quantiles([x[1] for x in cs], [0.5])[0]
                n_med = quantiles([x[0] for x in ns], [0.5])[0]
                n_medz = quantiles([x[1] for x in ns], [0.5])[0]
                sens = sum(1 for e, z, _n in cs if passes(e, z)) / float(len(cs))
                hits = [x for x in ns if passes(x[0], x[1])]
                A = auc([x[0] for x in cs], [x[0] for x in ns])
                print("%-7d %-8d | %7.2f %7.2f %6.0f%% | %7.2f %7.2f | %5.3f %8d %6d"
                      % (inner, width, c_med, c_medz, 100 * sens, n_med, n_medz,
                         A, len(hits), len(cs)))
                # rank by SENSITIVITY at the fixed thresholds, tie-broken by fewer null
                # pairs clearing them. Never by a threshold fitted to the data.
                key = (sens, -len(hits))
                if best is None or key > best[0]:
                    best = (key, inner, width, hits, sens, A)
        if best is None:
            print("\nNothing scoreable at any swept window.")
            return
        (bsens, _), bi, bw, bhits, sens, A = best
        print("\nBest sensitivity: inner=%d width=%d -- %.0f%% of KNOWN junctions clear "
              "enr>=%.1f and z>=%.1f (AUC %.3f)" % (bi, bw, 100 * bsens,
                                                    a.min_enrichment, a.min_z, A))
        if bsens < 0.5:
            print("\nFewer than half the known junctions clear the thresholds at ANY window")
            print("swept here, so this library cannot reliably resolve junctions. Report")
            print("that. Do NOT lower --min-z or --min-enrichment until something passes:")
            print("the controls are the calibration, and they say the test is blind.")
            return
        print("Inter-scaffold pairs clearing the same thresholds -- the CANDIDATES:")
        print("  %-36s %9s %9s" % ("pair", "enrich", "z"))
        for e, z, nm in sorted(bhits, key=lambda x: -x[1])[:20]:
            print("  %-36s %9.2f %9.2f" % (nm, e, z))
        if not bhits:
            print("  (none -- with %.0f%% sensitivity demonstrated, that is a real negative:"
                  % (100 * bsens))
            print("   Hi-C does not support any of these joins)")
        print("\nRe-run without --scan using --window-inner %d --window-bp %d to write "
              "calls." % (bi, bw))
        return

    ctrl_rows = [(c, r) for c, r in
                 ((c, score_control(c, a.window_inner, a.window_bp)) for c in ctrl) if r]
    ctrl_enr = [r["enr"] for _c, r in ctrl_rows if r["wcov"] >= a.min_window_cov]
    cmed, c10 = quantiles(ctrl_enr, [0.5, 0.1])

    rows = [(ka, kb, r) for ka, kb, r in
            ((ka, kb, score_pair(ka, kb, a.window_inner, a.window_bp))
             for i, ka in enumerate(cand) for kb in cand[i + 1:]) if r]
    null_enr = [r["enr"] for _a, _b, r in rows
                if r["tot"] >= a.min_total_contacts and r["wcov"] >= a.min_window_cov]
    nmed, n90 = quantiles(null_enr, [0.5, 0.9])
    # Thresholds are FIXED (--min-enrichment, --min-z). The controls do not set them --
    # they measure whether known junctions clear them, i.e. sensitivity. Deriving a
    # threshold from either distribution is unsound: the null contains the junctions being
    # searched for, and the controls are a stronger positive than a real missed junction.
    call_floor = a.min_enrichment
    ctrl_pass = [1 for _c, r in ctrl_rows
                 if r["wcov"] >= a.min_window_cov
                 and r["enr"] >= a.min_enrichment and r["z"] >= a.min_z]
    sensitivity = (len(ctrl_pass) / float(len(ctrl_enr))) if ctrl_enr else 0.0
    sensitive = sensitivity >= 0.5

    n_join = n_low = 0
    with open(a.out, "w") as fh:
        fh.write("# assembly\t%s\n" % a.assembly_id)
        fh.write("# params\tmin_scaffold_bp=%d\twindow_inner=%d\twindow_bp=%d\t"
                 "min_enrichment=%s\tmin_z=%s\tmin_total_contacts=%d\tmin_window_cov=%d\n"
                 % (a.min_scaffold_bp, a.window_inner, a.window_bp, a.min_enrichment,
                    a.min_z, a.min_total_contacts, a.min_window_cov))
        fh.write("# null\tinter_scaffold_pairs\t%d\tmedian\t%.3f\tp90\t%.3f\n"
                 % (len(null_enr), nmed, n90))
        fh.write("# positive_control\tartificial_splits\t%d\tmedian\t%.3f\tp10\t%.3f\n"
                 % (len(ctrl_enr), cmed, c10))
        fh.write("# sensitivity\t%s\t%.0f%%_of_known_junctions_clear_enr>=%.1f_and_z>=%.1f"
                 "\tnull_median\t%.3f\n"
                 % ("SENSITIVE" if sensitive else "INSENSITIVE", 100 * sensitivity,
                    a.min_enrichment, a.min_z, nmed))
        if not sensitive:
            fh.write("# WARNING\tfewer than half of the artificial splits -- KNOWN single "
                     "chromosomes -- clear these thresholds, so the test is blind at these "
                     "settings. All calls are suppressed to 'insensitive'. Do NOT lower "
                     "--min-z or --min-enrichment to make calls appear; the controls are "
                     "the calibration and they say the test cannot see.\n")
        for c, r in sorted(ctrl_rows, key=lambda x: -x[1]["enr"]):
            fh.write("# control\t%s\tenrichment\t%.3f\tz\t%.2f\tobs\t%d\texp\t%.1f\n"
                     % (c, r["enr"], r["z"], r["obs"], r["exp"]))
        fh.write("assembly\tscaf_a\tscaf_b\tlen_a\tlen_b\ttotal_contacts\tbest_corner\t"
                 "corner_contacts\texpected_contacts\tcorner_enrichment\tpoisson_z\t"
                 "window_coverage\tcontrol_median\tnull_p90\tcall\torientation\n")
        for ka, kb, r in sorted(rows, key=lambda x: -x[2]["enr"]):
            if not sensitive:
                call = "insensitive"
            elif r["tot"] < a.min_total_contacts or r["wcov"] < a.min_window_cov:
                call = "low_coverage"
                n_low += 1
            elif r["enr"] >= call_floor and r["z"] >= a.min_z:
                call = "join"
                n_join += 1
            else:
                call = "no_support"
            c = r["corner"]
            fh.write("\t".join(str(x) for x in (
                a.assembly_id, ka, kb, lens[ka], lens[kb], r["tot"], "%s-%s" % c,
                r["obs"], "%.2f" % r["exp"], "%.4f" % r["enr"], "%.2f" % r["z"],
                int(r["wcov"]), "%.3f" % cmed, "%.3f" % n90, call,
                "a_%s-b_%s" % c)) + "\n")

    sys.stderr.write(
        "[hic-junction %s] %d candidates, %d controls; sensitivity %.0f%% (control med "
        "%.2f, null med %.2f) -> %s; %d join, %d low_coverage\n"
        % (a.assembly_id, len(cand), len(ctrl_enr), 100 * sensitivity, cmed, nmed,
           "SENSITIVE" if sensitive else "INSENSITIVE", n_join, n_low))


if __name__ == "__main__":
    main()
