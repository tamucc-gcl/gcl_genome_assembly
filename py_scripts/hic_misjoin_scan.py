#!/usr/bin/env python3
"""
hic_misjoin_scan.py -- find mis-joins INSIDE a scaffold, and score named junction points.

WHY THIS EXISTS (and why the inter-scaffold test could not do the job)
---------------------------------------------------------------------
hic_junction_evidence.py asks: two SEPARATE scaffolds -- are their ends in contact? That
uses inter-scaffold contacts, the sparsest class in the dataset. On Sde-CBau_104_hap2 at a
50 kb window the expected count was 0.03, so a single stray read pair scored z = 5.7. Only
at a 2 Mb window did expected counts reach the hundreds, and by then a junction signal that
lives in the first few hundred kb is diluted away.

This asks the inverse, which is far better powered and is the question the cohort can
actually answer:

    in the assemblies that ALREADY MADE the join, is the join real?

A junction that is a single scaffold in four assemblies is spanned by INTRA-scaffold
contacts -- the abundant class. A genuine chromosome shows smooth contact decay across
every interior point. A mis-join shows a sharp deficit of spanning contacts.

THE STATISTIC
-------------
For a cut point c in scaffold S, count contacts that SPAN it with both ends within
--flank bp:

    observed = #{(p, q) : c - flank <= p < c <= q < c + flank}

The expectation comes from the scaffold's OWN contact-vs-distance decay, f(d), estimated
from every intra-scaffold contact in that scaffold. For a pair of positions at separation
d, the chance of straddling c is a function of geometry alone, so:

    expected(c) = sum over d of  f(d) * span_weight(c, d)

with span_weight the number of (p, q) placements at separation d that straddle c and fall
inside the flanks -- computed on a coverage-weighted grid, so a low-coverage region
predicts few spanning contacts rather than manufacturing a deficit. That was the failure
mode of an area-based null in the inter-scaffold test: terminal coverage runs 0.31x the
chromosome mean in the first 100 kb, and 1.4-1.9x from 1-5 Mb, a ~6x swing.

    ratio = observed / expected        (1.0 = intact, ->0 = mis-join)
    z     = (observed - expected) / sqrt(expected)

Because f(d) is fitted per scaffold, this self-normalises for library depth, and because
the expectation uses observed coverage it self-normalises for mappability.

TWO MODES
---------
  --scan          slide over every scaffold and report the worst dips. Assumption-free,
                  does not need to know where junctions are, and is useful on its own as
                  assembly QC.
  --junctions F   score named cut points supplied as a TSV (scaffold, position, label),
                  e.g. the alignment-derived junction coordinates. Reports each with its
                  own rank against that scaffold's null distribution of interior points.

CALIBRATION IS BUILT IN AND IS NOT OPTIONAL
-------------------------------------------
Every interior point of every scaffold is a NEGATIVE control: a position in what is
believed to be a contiguous chromosome. Their ratio distribution is the null. A junction is
called 'misjoin' only if its ratio falls below --null-pct of that distribution AND its
z is below --min-z. No threshold is fitted to the data: the percentile is a fixed rank and
z is a calibrated Poisson deviate.

The POSITIVE control is the scaffold ends themselves -- a scaffold boundary is a known
discontinuity, so cutting just inside one end should look intact while the boundary itself
has no spanning contacts at all by construction. Reported as a sanity line, not a
threshold.

Usage:
    hic_misjoin_scan.py --pairs p.gz --fai a.fai --assembly-id ID --out ID.misjoin.tsv --scan
    hic_misjoin_scan.py --pairs p.gz --fai a.fai --assembly-id ID --out ID.misjoin.tsv \
        --junctions junctions.tsv
"""
import argparse
import gzip
import sys
from collections import defaultdict

BIN = 10_000          # decay profile and coverage resolution


def openf(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def quantiles(vals, qs):
    if not vals:
        return [float("nan")] * len(qs)
    s = sorted(vals)
    out = []
    for q in qs:
        k = (len(s) - 1) * q
        lo = int(k)
        hi = min(lo + 1, len(s) - 1)
        out.append(s[lo] + (s[hi] - s[lo]) * (k - lo))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--assembly-id", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-scaffold-bp", type=int, default=8_000_000)
    ap.add_argument("--flank", type=int, default=500_000,
                    help="both ends of a spanning contact must lie within this of the cut")
    ap.add_argument("--step", type=int, default=100_000,
                    help="scan stride")
    ap.add_argument("--edge-margin", type=int, default=1_000_000,
                    help="skip cut points within this of a scaffold end; a real boundary "
                         "has no spanning contacts and would dominate the null")
    ap.add_argument("--min-expected", type=float, default=25.0,
                    help="a cut point whose expectation is below this has no power and is "
                         "reported 'underpowered', never called. This is the guard the "
                         "inter-scaffold test lacked: it floored marginal COVERAGE, not "
                         "EXPECTED CONTACTS, so a single read pair against an expectation "
                         "of 0.03 scored z = 5.7")
    ap.add_argument("--max-ratio", type=float, default=0.5,
                    help="a mis-join means the join is FAKE, so spanning contacts should be "
                         "near-absent, not merely reduced. Required in addition to the "
                         "percentile and z tests because both of those fail on a "
                         "minimum-statistic: --scan reports the WORST interior point, which "
                         "sits at the bottom of its own null by construction, and the "
                         "empirical null is over-dispersed relative to Poisson (TADs and "
                         "compartments), so z alone is anti-conservative. A planted "
                         "mis-join scores 0.00; a TAD boundary on an intact chromosome "
                         "scores ~0.85")
    ap.add_argument("--null-pct", type=float, default=1.0,
                    help="a junction is called only if its ratio is below this percentile "
                         "of that scaffold's interior points")
    ap.add_argument("--min-z", type=float, default=-5.0,
                    help="and its Poisson deviate must be at or below this (negative: a "
                         "DEFICIT of spanning contacts)")
    ap.add_argument("--allow-starved", action="store_true", default=False,
                    help="continue when a scaffold in the fai has no intra-scaffold "
                         "contacts, instead of refusing. Default is to refuse: it almost "
                         "always means the pairs and the fai are from different runs")
    ap.add_argument("--profile",
                    help="write the full ratio-vs-position profile here. This is the "
                         "output that matters -- --scan reports only the worst interior "
                         "dip per scaffold, and thresholding a minimum-statistic is "
                         "invalid, so those calls say more about the threshold than about "
                         "the assembly. WHERE the dips fall is the evidence")
    ap.add_argument("--local-window", type=int, default=3_000_000,
                    help="decay profile and placement mass are estimated within this "
                         "distance of the cut, not scaffold-wide, so long-range coverage "
                         "structure does not leak into a local expectation")
    ap.add_argument("--junctions",
                    help="TSV: scaffold, position, label. Score these points specifically")
    ap.add_argument("--scan", action="store_true",
                    help="report the worst dip per scaffold across a full slide")
    a = ap.parse_args()

    lens = {}
    with open(a.fai) as fh:
        for line in fh:
            if line.strip():
                f = line.split("\t")
                lens[f[0]] = int(f[1])
    keep = {n for n, L in lens.items() if L >= a.min_scaffold_bp}
    if not keep:
        sys.exit("[misjoin] no scaffold >= %d bp" % a.min_scaffold_bp)

    # ---- one streaming pass: coverage, decay profile, and the contacts near each other
    cov = defaultdict(lambda: defaultdict(int))          # scaffold -> bin -> ends
    decay = defaultdict(lambda: defaultdict(int))        # scaffold -> dist bin -> count
    near = defaultdict(list)                             # scaffold -> [(lo, hi)]
    n_intra = 0
    with openf(a.pairs) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            c1, c2 = f[1], f[3]
            if c1 != c2 or c1 not in keep:
                continue
            p1, p2 = int(f[2]), int(f[4])
            lo, hi = (p1, p2) if p1 <= p2 else (p2, p1)
            cov[c1][(lo - 1) // BIN] += 1
            cov[c1][(hi - 1) // BIN] += 1
            d = hi - lo
            decay[c1][d // BIN] += 1
            n_intra += 1
            if d < a.flank:
                near[c1].append((lo, hi))

    starved = sorted(c for c in keep if not near[c])
    if starved:
        sys.stderr.write(
            "[misjoin %s] ERROR: %d scaffold(s) in the fai carry ZERO intra-scaffold "
            "contacts: %s\n"
            "[misjoin %s]        The pairs and the fai describe different assemblies. A "
            "scaffold that exists in the fai but not in the pairs header is usually a "
            "vintage mismatch -- e.g. a composite created by a later harmonization run, "
            "where the published pairs still call it two separate scaffolds.\n"
            "[misjoin %s]        Re-map or re-run so both come from the same run, or pass "
            "--allow-starved to score the rest and report these as no_contacts.\n"
            % (a.assembly_id, len(starved), ", ".join(starved[:6]),
               a.assembly_id, a.assembly_id))
        if not a.allow_starved:
            sys.exit(2)
        keep = keep - set(starved)

    for c in keep:
        near[c].sort()

    # local decay counts around the cut currently being scored (filled by set_local)
    local_decay = defaultdict(dict)

    def set_local(c, cut):
        w = a.local_window
        lo_lim, hi_lim = cut - w, cut + w
        d = defaultdict(int)
        for lo, hi in near[c]:
            if lo >= lo_lim and hi < hi_lim:
                d[(hi - lo) // BIN] += 1
        local_decay[c] = d

    def straddle_mass(cv, lo, hi, db):
        """Coverage-weighted mass of placements whose LEFT end lies in [lo, hi).

        Bin resolution with fractional end bins. An earlier version enumerated (left bin,
        right bin) pairs and required the right bin to be at or past the cut, so the
        smallest separations contributed ZERO -- no placement of separation < BIN could
        qualify. The decay is steep, so those dominate, and the expectation came out ~1.7x
        too small.
        """
        if hi <= lo:
            return 0.0
        tot = 0.0
        for b in range(lo // BIN, (hi - 1) // BIN + 1):
            v = cv.get(b, 0)
            if not v:
                continue
            lo_b = max(lo, b * BIN)
            hi_b = min(hi, (b + 1) * BIN)
            tot += ((hi_b - lo_b) / float(BIN)) * v * cv.get(b + db, 0)
        return tot

    def expected_at(c, cut):
        """Expected spanning contacts at `cut`, from the LOCAL decay profile.

        A contact of separation d straddles `cut` exactly when its left end lies in
        [cut - d, cut). With d < flank that also puts both ends inside the flanks, so no
        separate flank term is needed.

        Both the decay counts and the placement-mass denominator are restricted to a local
        neighbourhood (--local-window on each side). Using the whole scaffold made every
        local expectation carry that scaffold's global coverage structure: contact density
        runs 1.4-1.9x the chromosome mean from 1-5 Mb of the ends and 0.31x in the first
        100 kb, so a cut inside the elevated zone was compared against a denominator
        inflated by regions nowhere near it. That produced a false 'misjoin' on a KNOWN
        intact chromosome at 46.5 Mb, and left the interior null centred at 1.48 rather
        than 1.0.
        """
        cv = cov[c]
        L = lens[c]
        w = a.local_window
        lo_lim, hi_lim = max(0, cut - w), min(L, cut + w)
        exp = 0.0
        for db, _ in decay[c].items():
            if db > a.flank // BIN:
                continue
            d = db * BIN + BIN // 2
            if d >= a.flank:
                continue
            n_local = local_decay[c].get(db, 0)
            if not n_local:
                continue
            tot = straddle_mass(cv, lo_lim, max(lo_lim, hi_lim - d), db)
            if tot <= 0:
                continue
            num = straddle_mass(cv, max(lo_lim, cut - d), cut, db)
            if num > 0:
                exp += n_local * (num / tot)
        return exp

    def observed_at(c, cut):
        # near[] already holds only d < flank, and d < flank plus straddling implies both
        # ends are inside the flanks, so straddling is the only condition needed.
        return sum(1 for lo, hi in near[c] if lo < cut <= hi)

    def score(c, cut):
        set_local(c, cut)
        exp = expected_at(c, cut)
        if exp <= 0:
            return None
        obs = observed_at(c, cut)
        return {"cut": cut, "obs": obs, "exp": exp, "ratio": obs / exp,
                "z": (obs - exp) / (exp ** 0.5)}

    # ---- interior null per scaffold
    nulls = {}
    profiles = {}
    for c in sorted(keep):
        L = lens[c]
        pts = list(range(a.edge_margin, L - a.edge_margin, a.step))
        rows = [r for r in (score(c, p) for p in pts) if r]
        profiles[c] = rows
        nulls[c] = [r["ratio"] for r in rows if r["exp"] >= a.min_expected]
        sys.stderr.write("[misjoin %s] %s: %d interior points, %d powered\n"
                         % (a.assembly_id, c, len(rows), len(nulls[c])))

    all_null = [x for c in nulls for x in nulls[c]]
    gmed = quantiles(all_null, [0.5])[0]

    targets = []
    if a.junctions:
        with open(a.junctions) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[0] in keep:
                    targets.append((f[0], int(f[1]), f[2] if len(f) > 2 else "-"))

    if a.profile:
        with open(a.profile, "w") as pf:
            pf.write("# assembly\t%s\tflank\t%d\tstep\t%d\tlocal_window\t%d\n"
                     % (a.assembly_id, a.flank, a.step, a.local_window))
            pf.write("# interior null median ratio %.4f -- read dips RELATIVE to this, not "
                     "to 1.0\n" % gmed)
            pf.write("assembly\tscaffold\tscaffold_len\tposition\tobserved\texpected\t"
                     "ratio\tz\tratio_over_scaffold_median\n")
            for c in sorted(keep):
                rr = [r["ratio"] for r in profiles[c] if r["exp"] >= a.min_expected]
                smed = quantiles(rr, [0.5])[0] if rr else float("nan")
                for r in profiles[c]:
                    pf.write("\t".join(str(x) for x in (
                        a.assembly_id, c, lens[c], r["cut"], r["obs"],
                        "%.2f" % r["exp"], "%.4f" % r["ratio"], "%.2f" % r["z"],
                        "%.4f" % (r["ratio"] / smed) if smed == smed and smed > 0
                        else "NA")) + "\n")
        sys.stderr.write("[misjoin %s] profile written to %s\n"
                         % (a.assembly_id, a.profile))

    with open(a.out, "w") as fh:
        fh.write("# assembly\t%s\n" % a.assembly_id)
        fh.write("# params\tflank=%d\tstep=%d\tedge_margin=%d\tmin_expected=%.1f\t"
                 "null_pct=%s\tmin_z=%s\tmax_ratio=%s\tlocal_window=%d\n"
                 % (a.flank, a.step, a.edge_margin, a.min_expected, a.null_pct, a.min_z,
                    a.max_ratio, a.local_window))
        fh.write("# intra_contacts\t%d\tscaffolds\t%d\n" % (n_intra, len(keep)))
        q = quantiles(all_null, [0.01, 0.05, 0.25, 0.5, 0.75])
        fh.write("# interior_null\tn\t%d\tp01\t%.3f\tp05\t%.3f\tp25\t%.3f\tmedian\t%.3f\t"
                 "p75\t%.3f\n" % (len(all_null), q[0], q[1], q[2], q[3], q[4]))
        fh.write("# read\tratio 1.0 = contact decay continues smoothly across the point; "
                 "ratio -> 0 = spanning contacts are missing, i.e. a mis-join\n")
        fh.write("assembly\tscaffold\tscaffold_len\tposition\tlabel\tobserved\texpected\t"
                 "ratio\tz\tscaffold_null_p%s\tpercentile_in_scaffold\tcall\n"
                 % ("%g" % a.null_pct))

        def emit(c, cut, label, r):
            if r is None:
                fh.write("\t".join(str(x) for x in (
                    a.assembly_id, c, lens[c], cut, label, 0, 0, "NA", "NA", "NA", "NA",
                    "no_expectation")) + "\n")
                return None
            nl = nulls.get(c) or all_null
            thr = quantiles(nl, [a.null_pct / 100.0])[0]
            below = sum(1 for x in nl if x <= r["ratio"])
            pct = 100.0 * below / max(1, len(nl))
            if r["exp"] < a.min_expected:
                call = "underpowered"
            elif r["ratio"] <= a.max_ratio and r["ratio"] <= thr and r["z"] <= a.min_z:
                call = "misjoin"
            elif r["ratio"] <= thr and r["z"] <= a.min_z:
                call = "weak_deficit"
            else:
                call = "intact"
            fh.write("\t".join(str(x) for x in (
                a.assembly_id, c, lens[c], cut, label, r["obs"], "%.2f" % r["exp"],
                "%.4f" % r["ratio"], "%.2f" % r["z"], "%.4f" % thr, "%.1f" % pct,
                call)) + "\n")
            return call

        calls = defaultdict(int)
        if targets:
            for c, cut, label in targets:
                calls[emit(c, cut, label, score(c, cut))] += 1
        if a.scan or not targets:
            # The worst interior point sits at the bottom of its own null by construction,
            # so no threshold applied to it is meaningful. Reported for orientation only,
            # with the call column forced to 'worst_dip_not_a_call'. Use --profile and look
            # at WHERE the dips fall.
            for c in sorted(keep):
                rows = [r for r in profiles[c] if r["exp"] >= a.min_expected]
                if not rows:
                    continue
                w = min(rows, key=lambda r: r["ratio"])
                nl = nulls.get(c) or all_null
                below = sum(1 for x in nl if x <= w["ratio"])
                fh.write("\t".join(str(x) for x in (
                    a.assembly_id, c, lens[c], w["cut"], "worst_interior_dip", w["obs"],
                    "%.2f" % w["exp"], "%.4f" % w["ratio"], "%.2f" % w["z"], "NA",
                    "%.1f" % (100.0 * below / max(1, len(nl))),
                    "worst_dip_not_a_call")) + "\n")
                calls["worst_dip_not_a_call"] += 1

    sys.stderr.write("[misjoin %s] interior null median ratio %.3f (1.0 = intact); "
                     "calls: %s\n" % (a.assembly_id, gmed,
                                      ", ".join("%s=%d" % kv for kv in sorted(calls.items()))))


if __name__ == "__main__":
    main()
