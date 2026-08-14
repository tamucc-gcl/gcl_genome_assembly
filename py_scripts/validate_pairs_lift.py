#!/usr/bin/env python3
"""
validate_pairs_lift.py -- prove the harmonized pairs lift is exact before anything depends
on it.

TEST A -- ARITHMETIC (exact, must pass)
    Compares the original and lifted pairs files, which share a lineage, so any difference
    is the lift's own doing. No aligner involved, so every check is exact:

      A1  pair count preserved (minus explicitly dropped)
      A2  per-scaffold end counts preserved under the rename
      A3  EXACT position moments -- n, sum(pos), sum(pos^2), min, max -- against the closed
          form for p' = L - p + 1. Base-resolution and O(1) memory.
      A4  binned profile shape (secondary; catches local reordering)
      A5  strand composition preserved for fwd, exchanged (+ <-> -) for rev
      A6  every lifted pair satisfies the upper-triangular convention in the new order
      A7  every lifted position within [1, length]

    A3 is the load-bearing check. A binned profile CANNOT detect a 1 bp shift -- it only
    moves a handful of reads across bin boundaries, which is indistinguishable from noise.
    sum(pos) changes by exactly n per bp of shift, so it cannot hide. This was found the
    hard way: an injected off-by-one produced a bin diff of 2 out of 200,000 ends.

TEST B -- CROSS-PATH (tolerance, informational)
    Compares a .cool built from the lifted pairs against the existing directly-mapped .cool
    for the same assembly. Both are in finalized coordinates, so they should agree -- but
    NOT bitwise: the two paths align to different (revcomp'd, renamed, reordered) references,
    so multi-mapping placement and MAPQ ties can differ slightly. Reports per-block relative
    differences and lets you judge.

    Requires cooler. Skipped if --cool-lifted / --cool-direct are not both given.

Usage:
    # Test A only
    validate_pairs_lift.py --orig-pairs pre.pairs.gz --lifted-pairs post.pairs.gz \
        --name-map nm.tsv --bin 100000

    # Test A + B
    validate_pairs_lift.py --orig-pairs pre.pairs.gz --lifted-pairs post.pairs.gz \
        --name-map nm.tsv --cool-lifted lifted.cool --cool-direct direct.cool
"""
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


def openf(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_name_map(path):
    lift = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        I = {k: hdr.index(k) for k in
             ("old_name", "new_name", "orient", "order", "length") if k in hdr}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            lift[f[I["old_name"]]] = (f[I["new_name"]],
                                      f[I["orient"]].strip().lower() == "rev",
                                      int(f[I["length"]]),
                                      int(f[I["order"]]))
    return lift


def identity_map(pairs_path):
    """Build an identity lift from a pairs header, for the round-trip check.

    A round trip (invert then forward) must return the input, so the expected relation
    between the two files is the identity: same names, orient fwd, same lengths. Reading
    the header means no map file is needed and nothing can be mis-specified.
    """
    lift = {}
    with openf(pairs_path) as fh:
        for i, line in enumerate(fh):
            if not line.startswith("#"):
                break
            if line.lower().startswith("#chromsize:"):
                p = line.split(":", 1)[1].split()
                if len(p) >= 2:
                    lift[p[0]] = (p[0], False, int(p[1]), len(lift) + 1)
    return lift


def read_chromsizes(path):
    """(name, length) pairs from a pairs file's #chromsize: header. Cheap -- header only."""
    cs = {}
    with openf(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.lower().startswith("#chromsize:"):
                p = line.split(":", 1)[1].split()
                if len(p) >= 2:
                    cs[p[0]] = int(p[1])
    return cs


def run_self_roundtrip(orig_pairs, name_map, workdir):
    """Run invert-then-forward internally, into a private temp dir.

    Exists because a paste-able command sequence without `set -e` leaves a stale
    intermediate behind when a lift REFUSES, and the next command silently reads it. That
    failure has happened twice and produced pages of misleading output both times. Running
    both passes here makes it impossible: unique paths, and a refusal aborts.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    lifter = os.path.join(here, "lift_harmonized_pairs.py")
    if not os.path.exists(lifter):
        sys.exit(f"--self-roundtrip needs lift_harmonized_pairs.py beside this script "
                 f"({lifter} not found)")
    back = os.path.join(workdir, "back.pairs")
    fwd = os.path.join(workdir, "roundtrip.pairs")
    for label, args in (("invert", ["--pairs", orig_pairs, "--invert", "--out", back]),
                        ("forward", ["--pairs", back, "--out", fwd])):
        cmd = [sys.executable, lifter, "--name-map", name_map,
               "--chrom-sizes", os.path.join(workdir, f"{label}.chrom.sizes"),
               "--stats", os.path.join(workdir, f"{label}.stats.tsv")] + args
        sys.stderr.write(f"[self-roundtrip] {label} pass\n")
        r = subprocess.run(cmd)
        if r.returncode != 0:
            sys.exit(f"[self-roundtrip] the {label} pass refused (exit {r.returncode}). "
                     f"Nothing was compared. The pairs and the map are from different "
                     f"assemblies -- see the message above.")
    return fwd


def scan(path, binsize, order=None, rev_len=None):
    """Single streaming pass. O(scaffolds) memory -- never stores per-pair data.

    rev_len maps a reversed scaffold to its length; those scaffolds are binned from the FAR
    END, ((L - p) // B), so the profile matches the original's ((p - 1) // B) exactly for
    any length. Binning both from the start only lines up when the bin size divides the
    length, which on real assemblies is never -- the check reported 0/0 and called it a
    pass, which is worse than not running it.

    Per scaffold: end count, sum(pos), sum(pos^2), min, max, strand counts, binned profile.
    The moments are the load-bearing checks: they are exact at base resolution and O(1)
    memory, so a uniform 1 bp shift changes sum(pos) by exactly n and cannot hide. A binned
    profile cannot do this -- a 1 bp shift only moves reads across bin boundaries, which is
    a handful of reads out of millions.
    """
    st = {}
    n = 0
    tri_bad = 0

    def slot(c):
        v = st.get(c)
        if v is None:
            v = st[c] = {"n": 0, "s": 0, "s2": 0, "min": None, "max": 0,
                         "plus": 0, "minus": 0, "prof": defaultdict(int)}
        return v

    with openf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            n += 1
            c1, p1, c2, p2 = f[1], int(f[2]), f[3], int(f[4])
            for c, p, sd in ((c1, p1, f[5]), (c2, p2, f[6])):
                v = slot(c)
                v["n"] += 1
                v["s"] += p
                v["s2"] += p * p
                if v["min"] is None or p < v["min"]:
                    v["min"] = p
                if p > v["max"]:
                    v["max"] = p
                if sd == "+":
                    v["plus"] += 1
                elif sd == "-":
                    v["minus"] += 1
                Lr = rev_len.get(c) if rev_len else None
                v["prof"][((Lr - p) if Lr is not None else (p - 1)) // binsize] += 1
            if order is not None:
                i1 = order.get(c1, 0)
                i2 = order.get(c2, 0)
                if i1 > i2 or (i1 == i2 and p1 > p2):
                    tri_bad += 1
    return n, st, tri_bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orig-pairs", required=True)
    ap.add_argument("--lifted-pairs",
                    required=not any(x == "--self-roundtrip" for x in sys.argv),
                    help="not needed with --self-roundtrip, which produces it internally")
    ap.add_argument("--name-map",
                    help="applied_lift map relating orig -> lifted. Not needed with "
                         "--roundtrip, where the names are unchanged by construction")
    ap.add_argument("--self-roundtrip", action="store_true", default=False,
                    help="run the whole round trip internally: invert then forward into a "
                         "private temp dir, then validate. Needs --orig-pairs and "
                         "--name-map only. Preferred over --roundtrip, because there is no "
                         "intermediate file that can go stale between commands.")
    ap.add_argument("--roundtrip", action="store_true", default=False,
                    help="orig and lifted are the two ends of an invert-then-forward round "
                         "trip, so the expected relation is the identity. Runs the same "
                         "exact checks with orient=fwd throughout. Uses only published "
                         "files -- no pre-finalize pairs required.")
    ap.add_argument("--bin", type=int, default=100_000,
                    help="bin size for the marginal profile test; must divide evenly into "
                         "each scaffold length for the reversal check to be exact, so the "
                         "test reports exact-bin scaffolds separately")
    ap.add_argument("--cool-lifted")
    ap.add_argument("--cool-direct")
    ap.add_argument("--cool-tolerance", type=float, default=0.02,
                    help="max acceptable relative difference in per-block contact sums")
    a = ap.parse_args()

    fails, warns = [], []
    tmpdir = None
    if a.self_roundtrip:
        if not a.name_map:
            sys.exit("--self-roundtrip needs --name-map")
        tmpdir = tempfile.mkdtemp(prefix="rtcheck.")
        a.lifted_pairs = run_self_roundtrip(a.orig_pairs, a.name_map, tmpdir)
        a.roundtrip = True
    if a.roundtrip:
        lift = identity_map(a.orig_pairs)
        if not lift:
            sys.exit("--roundtrip needs #chromsize: lines in the original pairs header")
        print(f"[roundtrip] identity map over {len(lift)} scaffolds from the pairs header")
    elif a.name_map:
        lift = read_name_map(a.name_map)
    else:
        sys.exit("give --name-map, or --roundtrip to compare an invert-then-forward pair")

    # ---- A0: compare the two files' #chromsize headers BEFORE scanning. Names alone are a
    # weak discriminator -- two assemblies of the same species share a naming scheme, so a
    # stale file from a sibling assembly scored 97% on a name-only check. Lengths settle it,
    # and reading headers costs milliseconds instead of scanning 180M lines to reach a
    # conclusion that was available up front. Mismatched inputs are an INVALID TEST, not a
    # failed lift, so this exits 2 and runs nothing else.
    want_cs = {new: L for _o, (new, _r, L, _rk) in lift.items()}
    got_cs = read_chromsizes(a.lifted_pairs)
    if got_cs and want_cs != got_cs:
        only_want = sorted(set(want_cs) - set(got_cs))
        only_got = sorted(set(got_cs) - set(want_cs))
        badlen = sorted(c for c in set(want_cs) & set(got_cs) if want_cs[c] != got_cs[c])
        print("=" * 78)
        print("MISMATCHED INPUTS -- nothing was compared")
        print("=" * 78)
        print(f"  expected {len(want_cs)} scaffolds, lifted file has {len(got_cs)}")
        if only_want:
            print(f"  absent from lifted   ({len(only_want)}): {only_want[:4]}")
        if only_got:
            print(f"  unexpected in lifted ({len(only_got)}): {only_got[:4]}")
        if badlen:
            print(f"  length disagreements ({len(badlen)}), e.g. "
                  + "; ".join(f"{c} expected {want_cs[c]:,} got {got_cs[c]:,}"
                              for c in badlen[:3]))
        print()
        print("  Same-species assemblies share a naming scheme, so matching NAMES prove")
        print("  nothing -- the lengths are what settle it. Most likely an earlier lift")
        print("  REFUSED and a leftover file from a previous command was read instead.")
        print("  Use --self-roundtrip, which runs both passes into a private temp dir and")
        print("  cannot pick up a stale intermediate.")
        return 2

    order = {new: o for _old, (new, _r, _L, o) in lift.items()}
    sys.stderr.write("[validate] scanning original pairs...\n")
    n0, st0, _ = scan(a.orig_pairs, a.bin)
    sys.stderr.write("[validate] scanning lifted pairs...\n")
    rev_len = {new: L for _o, (new, rv, L, _r) in lift.items() if rv}
    n1, st1, tri_bad = scan(a.lifted_pairs, a.bin, order, rev_len)

    print("=" * 78)
    print("TEST A -- ARITHMETIC (exact)")
    print("=" * 78)

    ok = n0 == n1
    print(f"A1 pair count            orig={n0:,}  lifted={n1:,}  {'PASS' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"A1 pair count {n0} -> {n1}")

    bad2 = []
    for old, (new, rev, L, _o) in lift.items():
        e0 = st0.get(old, {}).get("n", 0)
        e1 = st1.get(new, {}).get("n", 0)
        if e0 != e1:
            bad2.append((old, new, e0, e1))
    print(f"A2 per-scaffold ends     {len(lift) - len(bad2)}/{len(lift)} match  "
          f"{'PASS' if not bad2 else 'FAIL'}")
    for b in bad2[:5]:
        print(f"     {b[0]} -> {b[1]}: {b[2]:,} vs {b[3]:,}")
    if bad2:
        fails.append(f"A2 {len(bad2)} scaffolds with changed end counts")

    # A3 -- exact position moments. For orient=rev the map is p' = L - p + 1, so:
    #        sum(p')   = n(L+1) - sum(p)
    #        sum(p'^2) = n(L+1)^2 - 2(L+1)sum(p) + sum(p^2)
    #        min(p')   = L - max(p) + 1        max(p') = L - min(p) + 1
    #      A uniform off-by-one changes sum(p') by exactly n. Nothing hides here.
    bad3 = []
    checked3 = 0
    for old, (new, rev, L, _o) in lift.items():
        v0, v1 = st0.get(old), st1.get(new)
        if not v0 or not v0["n"]:
            continue
        checked3 += 1
        if v1 is None:
            bad3.append((old, new, "rev" if rev else "fwd", "absent from the lifted file"))
            continue
        n_e = v0["n"]
        if rev:
            K = L + 1
            want_s = n_e * K - v0["s"]
            want_s2 = n_e * K * K - 2 * K * v0["s"] + v0["s2"]
            want_min = L - v0["max"] + 1
            want_max = L - v0["min"] + 1
        else:
            want_s, want_s2 = v0["s"], v0["s2"]
            want_min, want_max = v0["min"], v0["max"]
        why = []
        if v1["s"] != want_s:
            d = v1["s"] - want_s
            shift = d / n_e if n_e else 0
            why.append("sum(pos) off by %+d (= %+.3f bp/read)" % (d, shift))
        if v1["s2"] != want_s2:
            why.append("sum(pos^2) mismatch")
        if v1["min"] != want_min or v1["max"] != want_max:
            why.append("min/max %d..%d expected %d..%d"
                       % (v1["min"], v1["max"], want_min, want_max))
        if why:
            bad3.append((old, new, "rev" if rev else "fwd", "; ".join(why)))
    print(f"A3 position moments      {checked3 - len(bad3)}/{checked3} exact  "
          f"{'PASS' if not bad3 else 'FAIL'}")
    for b in bad3[:6]:
        print(f"     {b[0]} -> {b[1]} [{b[2]}]: {b[3]}")
    if bad3:
        fails.append(f"A3 {len(bad3)} scaffolds with wrong positions "
                     f"(a whole-bp shift means L - pos + 1 is wrong)")

    # A4 -- binned profile shape. Secondary to A3: catches local reordering that preserves
    # the moments. Reversed scaffolds are binned from the far end in scan(), so the two
    # profiles are directly comparable at any length.
    exact4, checked4, bad4 = 0, 0, []
    for old, (new, rev, L, _o) in lift.items():
        v0, v1 = st0.get(old), st1.get(new)
        if not v0 or not v0["n"]:
            continue
        checked4 += 1
        if v1 is None:
            bad4.append((old, new, "rev" if rev else "fwd", -1))
            continue
        if all(v0["prof"].get(b, 0) == v1["prof"].get(b, 0)
               for b in set(v0["prof"]) | set(v1["prof"])):
            exact4 += 1
        else:
            bad4.append((old, new, "rev" if rev else "fwd",
                         sum(abs(v0["prof"].get(b, 0) - v1["prof"].get(b, 0))
                             for b in set(v0["prof"]) | set(v1["prof"]))))
    verdict4 = "SKIP" if checked4 == 0 else ("PASS" if not bad4 else "FAIL")
    print(f"A4 profile shape         {exact4}/{checked4} exact  {verdict4}")
    for b in bad4[:5]:
        print(f"     {b[0]} -> {b[1]} [{b[2]}]: total abs bin diff {b[3]:,}")
    if bad4:
        fails.append(f"A4 {len(bad4)} scaffolds with wrong profile shape")
    if checked4 == 0:
        warns.append("A4 checked nothing")

    bad5 = []
    for old, (new, rev, L, _o) in lift.items():
        v0, v1 = st0.get(old), st1.get(new)
        if not v0:
            continue
        if v1 is None:
            bad5.append((old, new, "rev" if rev else "fwd", "absent", "-"))
            continue
        wp, wm = (v0["minus"], v0["plus"]) if rev else (v0["plus"], v0["minus"])
        if v1["plus"] != wp or v1["minus"] != wm:
            bad5.append((old, new, "rev" if rev else "fwd",
                         f"{v1['plus']}/{v1['minus']}", f"{wp}/{wm}"))
    print(f"A5 strand composition    {len(lift) - len(bad5)}/{len(lift)} match  "
          f"{'PASS' if not bad5 else 'FAIL'}")
    for b in bad5[:5]:
        print(f"     {b[0]} -> {b[1]} [{b[2]}]: got +/- {b[3]}, expected {b[4]}")
    if bad5:
        fails.append(f"A5 {len(bad5)} scaffolds with wrong strand composition")

    print(f"A6 upper triangle        {n1 - tri_bad:,}/{n1:,} conform  "
          f"{'PASS' if not tri_bad else 'FAIL'}")
    if tri_bad:
        fails.append(f"A6 {tri_bad} pairs violate the upper-triangular convention")

    lenof = {new: L for _old, (new, _r, L, _o) in lift.items()}
    bad7 = [(c, st1[c]["min"], st1[c]["max"], lenof.get(c)) for c in st1
            if lenof.get(c) is not None
            and (st1[c]["min"] < 1 or st1[c]["max"] > lenof[c])]
    print(f"A7 positions in range    {len(st1) - len(bad7)}/{len(st1)} scaffolds  "
          f"{'PASS' if not bad7 else 'FAIL'}")
    for b in bad7[:5]:
        print(f"     {b[0]}: {b[1]:,}..{b[2]:,} vs length {b[3]:,}")
    if bad7:
        fails.append(f"A7 {len(bad7)} scaffolds with out-of-range positions")

    # ---- TEST B
    if a.cool_lifted and a.cool_direct:
        print()
        print("=" * 78)
        print("TEST B -- CROSS-PATH (tolerance)")
        print("=" * 78)
        try:
            import cooler
            import numpy as np
        except ImportError:
            print("  SKIP: cooler/numpy unavailable in this environment")
        else:
            cl = cooler.Cooler(a.cool_lifted)
            cd = cooler.Cooler(a.cool_direct)
            sl, sd = dict(cl.chromsizes), dict(cd.chromsizes)
            same = sl == sd
            print(f"B1 chrom set + lengths   {len(sl)} vs {len(sd)}  "
                  f"{'PASS' if same else 'FAIL'}")
            if not same:
                only = set(sl) ^ set(sd)
                print(f"     symmetric difference: {sorted(only)[:6]}")
                fails.append("B1 chromsizes differ between the two cools")
            tl = cl.info.get("sum", 0)
            td = cd.info.get("sum", 0)
            rel = abs(tl - td) / max(1, td)
            print(f"B2 total contacts        lifted={tl:,} direct={td:,}  "
                  f"rel diff {rel:.4f}  "
                  f"{'PASS' if rel <= a.cool_tolerance else 'WARN'}")
            if rel > a.cool_tolerance:
                warns.append(f"B2 total contacts differ by {rel:.1%}")
            if same:
                chroms = list(sl)
                worst, rows = (0.0, None), []
                for i, ci in enumerate(chroms):
                    for cj in chroms[i:]:
                        vl = float(np.nansum(cl.matrix(balance=False,
                                                       sparse=True).fetch(ci, cj).data))
                        vd = float(np.nansum(cd.matrix(balance=False,
                                                       sparse=True).fetch(ci, cj).data))
                        r = abs(vl - vd) / max(1.0, vd)
                        rows.append((r, ci, cj, vl, vd))
                        if r > worst[0]:
                            worst = (r, (ci, cj))
                rows.sort(reverse=True)
                n_over = sum(1 for r, *_ in rows if r > a.cool_tolerance)
                print(f"B3 per-block sums        {len(rows) - n_over}/{len(rows)} within "
                      f"{a.cool_tolerance:.0%}  worst {worst[0]:.4f} at {worst[1]}  "
                      f"{'PASS' if not n_over else 'WARN'}")
                for r, ci, cj, vl, vd in rows[:6]:
                    print(f"     {ci} x {cj}: lifted={vl:,.0f} direct={vd:,.0f} rel={r:.4f}")
                if n_over:
                    warns.append(f"B3 {n_over} blocks over tolerance")

    if tmpdir:
        shutil.rmtree(tmpdir, ignore_errors=True)
    print()
    print("=" * 78)
    if fails:
        print("VERDICT: FAIL -- do not use the lift")
        for f in fails:
            print("  " + f)
        print("\nTest A failures are lift bugs, not aligner noise. A3 in particular means "
              "the coordinate reversal is wrong.")
        return 1
    print("VERDICT: Test A exact." + ("  Test B within tolerance." if a.cool_lifted
                                      and a.cool_direct and not warns else ""))
    for w in warns:
        print("  WARN: " + w)
    if warns:
        print("\nTest B differences are expected at some level -- the two paths align to "
              "different references, so multi-mapping placement and MAPQ ties can diverge. "
              "Judge whether the magnitude is acceptable; Test A is the correctness proof.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
