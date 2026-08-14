#!/usr/bin/env python3
"""
hic_junction_evidence.py -- does Hi-C say two scaffolds of ONE assembly are a single
chromosome that the scaffolder failed to join?

This is the only evidence channel that can break a tie the cross-assembly alignment
consensus cannot. It is a WITHIN-assembly measurement, so it is independent of the
harmonization frame, of the choice of reference, and of every other assembly in the cohort.

THE TEST
--------
If scaffolds A and B are two halves of one chromosome, Hi-C contacts between them
concentrate at the two ends that abut. If they are genuinely separate chromosomes, contacts
between them are flat inter-chromosomal background.

For each candidate scaffold pair, contacts are counted in the four terminal windows
(A-start/B-start, A-start/B-end, A-end/B-start, A-end/B-end) and compared to the density
over the whole A x B block:

    corner_enrichment = (corner_contacts / corner_area) / (total_AB / (len_A * len_B))

Self-normalizing: it asks whether contacts are CONCENTRATED at a corner, not how many there
are. A repeat-driven or mappability artifact elevates the whole block and cancels out; only
corner localization survives. Which of the four corners wins also gives the relative
ORIENTATION, which frame decoupling needs to reorient composites.

DECISION RULE, FIXED BEFORE LOOKING
-----------------------------------
Under the null that A and B are separate chromosomes, contacts are spread uniformly over the
A x B block, so the expected count in a corner is

    expected = total_AB * (window_a * window_b) / (len_a * len_b)

and the observed count is Poisson about it. A pair is called 'join' when its best corner
clears BOTH

    enrichment = observed / expected  >=  --min-corner-enrichment
    z          = (observed - expected) / sqrt(expected)  >=  --min-z

The z term is what makes this coverage-aware: a 2x enrichment on an expectation of 5 is
noise (z = 1.3), the same 2x on an expectation of 1000 is not (z = 31). Both terms are
needed -- z alone would call a 1.2x enrichment significant on a deep library, and
enrichment alone would call noise on a shallow one.

An earlier version used a percentile of the other pairs' enrichments as the threshold. That
was wrong: real junctions are a few percent of pairs, so they sit INSIDE the top 5% and the
p95 lands among the signal it is supposed to be measured against. On a synthetic test with
one planted junction out of six pairs the threshold rose to 7.8 and the junction at 10.0
only just cleared its own contribution; two junctions would have suppressed each other. The
Poisson null needs no background distribution at all. The median enrichment is still
reported, as a sanity check that the null sits near 1.

Pairs with fewer than --min-total-contacts are reported 'low_coverage' and never called
either way -- absence of contact is not evidence of separation.

WHY A DIRECT PASS AND NOT COOLER
--------------------------------
Bin-free. One streaming pass, O(candidate_pairs) counters, no resolution to choose and no
quantization of the window edges. It also removes a dependency and the need to pick a bin
size that is fine enough to localize a junction but coarse enough for stable counts.

CANDIDATE SCAFFOLDS
-------------------
Every scaffold >= --min-scaffold-bp in this assembly's own fai. Deliberately NOT the
drop-off caller used elsewhere: this must not depend on the harmonization frame, and
over-inclusion is harmless -- an extra scaffold contributes rows that get called
'no_support' or 'low_coverage' and widens the background distribution slightly.

Usage:
    hic_junction_evidence.py --pairs in.pairs.gz --fai asm.fa.fai \
        --assembly-id <id> --out <id>.junction_evidence.tsv
"""
import argparse
import gzip
import sys
from collections import defaultdict

CORNERS = (("start", "start"), ("start", "end"), ("end", "start"), ("end", "end"))


def openf(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def median(vals):
    if not vals:
        return 0.0
    s = sorted(vals)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True, help="pairs.gz in THIS assembly's coordinates")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--assembly-id", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-scaffold-bp", type=int, default=1_000_000,
                    help="candidate scaffolds are those at least this long")
    ap.add_argument("--window-frac", type=float, default=0.25,
                    help="terminal window as a fraction of each scaffold's length")
    ap.add_argument("--window-max-bp", type=int, default=10_000_000,
                    help="cap on the terminal window")
    ap.add_argument("--min-corner-enrichment", type=float, default=2.0,
                    help="floor on observed/expected corner contacts")
    ap.add_argument("--min-z", type=float, default=5.0,
                    help="floor on the Poisson deviate (obs - exp) / sqrt(exp). Guards "
                         "against calling noise on a shallow library; 5.0 over ~10^5 tests "
                         "is still conservative after multiple-testing correction")
    ap.add_argument("--min-total-contacts", type=int, default=100,
                    help="below this a pair is 'low_coverage' and never called either way")
    ap.add_argument("--max-candidates", type=int, default=1000,
                    help="safety valve: if an assembly has more candidates than this, keep "
                         "only the longest (a shattered assembly would otherwise generate "
                         "millions of pair counters for no useful signal)")
    a = ap.parse_args()

    # ---- candidates and their terminal windows
    lens = {}
    with open(a.fai) as fh:
        for line in fh:
            if line.strip():
                f = line.split("\t")
                lens[f[0]] = int(f[1])
    cand = sorted((n for n, L in lens.items() if L >= a.min_scaffold_bp),
                  key=lambda n: (-lens[n], n))
    trimmed = 0
    if len(cand) > a.max_candidates:
        trimmed = len(cand) - a.max_candidates
        cand = cand[:a.max_candidates]
    idx = {n: i for i, n in enumerate(cand)}
    win = {n: max(1, min(a.window_max_bp, int(a.window_frac * lens[n]))) for n in cand}

    sys.stderr.write("[hic-junction %s] %d candidate scaffolds (>= %d bp)%s\n" % (
        a.assembly_id, len(cand), a.min_scaffold_bp,
        "; %d dropped by --max-candidates" % trimmed if trimmed else ""))

    # ---- one streaming pass
    total = defaultdict(int)
    corner = defaultdict(lambda: defaultdict(int))
    n_pairs = n_inter = 0
    with openf(a.pairs) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_pairs += 1
            f = line.rstrip("\n").split("\t")
            c1, c2 = f[1], f[3]
            if c1 == c2:
                continue
            i1 = idx.get(c1)
            i2 = idx.get(c2)
            if i1 is None or i2 is None:
                continue
            n_inter += 1
            p1, p2 = int(f[2]), int(f[4])
            if i1 <= i2:
                ka, kb, pa, pb = c1, c2, p1, p2
            else:
                ka, kb, pa, pb = c2, c1, p2, p1
            key = (ka, kb)
            total[key] += 1
            ea = "start" if pa <= win[ka] else ("end" if pa > lens[ka] - win[ka] else None)
            if ea is None:
                continue
            eb = "start" if pb <= win[kb] else ("end" if pb > lens[kb] - win[kb] else None)
            if eb is None:
                continue
            corner[key][(ea, eb)] += 1

    # ---- enrichment per pair, then the per-assembly background
    rows = []
    for (ka, kb), tot in total.items():
        La, Lb = lens[ka], lens[kb]
        Wa, Wb = win[ka], win[kb]
        exp = tot * (float(Wa) * Wb) / (float(La) * Lb)
        best = (0.0, None, 0, 0.0)
        for c in CORNERS:
            cnt = corner[(ka, kb)].get(c, 0)
            enr = (cnt / exp) if exp > 0 else 0.0
            if enr > best[0]:
                z = (cnt - exp) / (exp ** 0.5) if exp > 0 else 0.0
                best = (enr, c, cnt, z)
        rows.append({"a": ka, "b": kb, "La": La, "Lb": Lb, "Wa": Wa, "Wb": Wb,
                     "tot": tot, "exp": exp, "enr": best[0], "corner": best[1],
                     "cnt": best[2], "z": best[3]})

    scored = [r for r in rows if r["tot"] >= a.min_total_contacts]
    med = median([r["enr"] for r in scored])

    n_join = n_low = 0
    with open(a.out, "w") as fh:
        fh.write("# assembly\t%s\n" % a.assembly_id)
        fh.write("# params\tmin_scaffold_bp=%d\twindow_frac=%s\twindow_max_bp=%d\t"
                 "min_corner_enrichment=%s\tmin_z=%s\tmin_total_contacts=%d\n"
                 % (a.min_scaffold_bp, a.window_frac, a.window_max_bp,
                    a.min_corner_enrichment, a.min_z, a.min_total_contacts))
        fh.write("# candidates\t%d\tpairs_scored\t%d\tmedian_enrichment\t%.4f"
                 "\t(null should sit near 1.0)\n" % (len(cand), len(scored), med))
        fh.write("# contacts\ttotal_pairs\t%d\tinter_candidate\t%d\n" % (n_pairs, n_inter))
        fh.write("assembly\tscaf_a\tscaf_b\tlen_a\tlen_b\twindow_a\twindow_b\t"
                 "total_contacts\tbest_corner\tcorner_contacts\texpected_contacts\t"
                 "corner_enrichment\tpoisson_z\tcall\torientation\n")
        for r in sorted(rows, key=lambda r: -r["z"]):
            if r["tot"] < a.min_total_contacts:
                call = "low_coverage"
                n_low += 1
            elif r["enr"] >= a.min_corner_enrichment and r["z"] >= a.min_z:
                call = "join"
                n_join += 1
            else:
                call = "no_support"
            c = r["corner"] or ("-", "-")
            fh.write("\t".join(str(x) for x in (
                a.assembly_id, r["a"], r["b"], r["La"], r["Lb"], r["Wa"], r["Wb"],
                r["tot"], "%s-%s" % c, r["cnt"], "%.2f" % r["exp"],
                "%.4f" % r["enr"], "%.2f" % r["z"],
                call, "a_%s-b_%s" % c if r["corner"] else "-")) + "\n")

    sys.stderr.write("[hic-junction %s] %d pairs scored, %d join, %d low_coverage; "
                     "median enrichment %.2f (null ~1.0); thresholds enr>=%.1f z>=%.1f\n"
                     % (a.assembly_id, len(scored), n_join, n_low, med,
                        a.min_corner_enrichment, a.min_z))


if __name__ == "__main__":
    main()
