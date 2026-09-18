#!/usr/bin/env python3
# ======================================================================================
# chimera_hic_pairs.py
#
# Contig-space Hi-C pairs -> pairs in a candidate scaffold's own coordinates, ready for
# `cooler cload pairs`.
#
# WHY THIS EXISTS RATHER THAN RE-MAPPING
# -------------------------------------
# Hi-C evidence needs a contact map in the coordinates of the assembly being cut, and at
# detection time nothing existing is in them:
#
#   FILTER_HIC_BAM             contig space -- no scaffold coordinates at all
#   FILTER_HIC_BAM_SCAFFOLD    round-1 scaffolds -- and round-1 `scaffold_1` is a DIFFERENT
#                              sequence from round-2 `scaffold_1`, same name
#   <asm>_final.mcool          produced by FINAL_HIC_MAPS, downstream of harmonization, so
#                              it does not exist yet on this run
#
# Relying on a PREVIOUS run's mcool was considered and rejected: on a fresh `chimera_break =
# auto` run there is no previous mcool, so Hi-C would be silently unavailable exactly when an
# unattended cut is being made.
#
# Re-mapping reads to a mini reference would work but needs bwa + samtools + pairtools +
# cooler chained together. Unnecessary: FILTER_HIC_BAM already publishes deduplicated
# unique-unique pairs in CONTIG space, and the AGP chain already says where every contig
# lands in the final scaffold. So this is a coordinate translation, not an alignment.
#
# THREE THINGS THE TRANSLATION MUST GET RIGHT
# -------------------------------------------
# 1. ORIENTATION. A contig placed '-' is reverse-complemented into its scaffold, so its
#    positions reflect: final = hi - (pos - 1) rather than lo + (pos - 1). And orientation
#    COMPOSES across the two AGP rounds -- a contig reversed into a round-1 scaffold that is
#    itself reversed into round-2 is forward overall.
#
# 2. UPPER TRIANGLE. The input header says `#shape: upper triangle`: pairs are stored with
#    (chrom1,pos1) <= (chrom2,pos2) in CONTIG order. Translation breaks that -- two contigs
#    adjacent in the scaffold can sort oppositely by name, and a reversed contig inverts
#    internally. Each pair has to be RE-CANONICALISED, swapping mates and strands together.
#    Getting this wrong yields a cool that loads but is silently asymmetric, and every
#    cross-contact ratio read from a half-empty triangle would be wrong in a plausible-
#    looking way.
#
# 3. SORT ORDER. The input is `#sorted: chr1-chr2-pos1-pos2` in contig order, which
#    translation destroys. `cooler cload pairs` needs sorted input, so the caller must sort
#    the output; this script does not, because sorting a few million lines belongs in the
#    shell and streaming here keeps memory flat over a ~1.9 GB gzipped input.
#
# WHAT IS CHECKED, AND WHAT IS NOT
# --------------------------------
# Two self-checks, and it is worth being precise about their reach because the first looks
# stronger than it is:
#
#   distance invariance   an intra-contig pair's distance must survive translation. This only
#                         proves the transform is RIGID. Reflecting from the wrong anchor
#                         (f_lo - off instead of f_hi - off) preserves every distance and a
#                         deliberately broken variant passed it 3/3.
#
#   range containment     every translated position must land inside the contig's own final
#                         interval. THIS is what catches a wrong anchor -- the broken variant
#                         produced -198 for an interval of [1, 1000] and aborts.
#
# NEITHER CATCHES A WRONG ORIENTATION. Using the round-1 orientation alone instead of
# composing it with round 2 still maps every position inside the correct interval, merely
# backwards within it, so both checks pass. That cannot be self-checked here: orientation is
# DERIVED from the two AGP levels, so there is nothing independent in the AGP to compare it
# against.
#
# It is instead verified two other ways. By hand on a fixture where sc_A = ctgF(-) + ctgE(+)
# is itself placed reversed, giving ctgF forward at 1101-2100 and ctgE reversed at 1-1000 --
# the two reversals cancelling for one contig and not the other. And by agp_joins.py's
# neighbour-consistency check, which exercises the same composition on 173 real joins of
# Sde-CTlk_104_hap2 scaffold_3 and caught a left/right inversion that reading the code had
# not.
#
# ONE PASS PER ASSEMBLY, not per candidate: the input is ~1.9 GB gzipped, so every candidate
# scaffold is filtered in the same stream.
#
# Dependencies: python3 stdlib.
#
# USAGE
#   chimera_hic_pairs.py --pairs contig.pairs.gz --round1 r1.agp [--round2 r2.agp] \
#       --scaffolds scaffold_1,scaffold_3 --outdir . --label <asm> [--min-mapq-pair UU]
# ======================================================================================

import argparse
import gzip
import os
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_agp(path):
    """AGP -> {object: [component dicts]}, sorted by part. Coordinates 1-based inclusive."""
    comps = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[4] in ("N", "U"):
                continue
            comps.setdefault(f[0], []).append(
                dict(part=int(f[3]), obeg=int(f[1]), oend=int(f[2]), cid=f[5],
                     cbeg=int(f[6]), cend=int(f[7]), orient=f[8]))
    for k in comps:
        comps[k].sort(key=lambda x: x["part"])
    return comps


def contig_placements(c1, c2, want):
    """contig -> (final_lo, final_hi, contig_sub_lo, orientation) for one final scaffold.

    Chains contig -> round-1 scaffold -> round-2 scaffold, composing orientation. A contig
    reversed into round 1 inside a round-1 scaffold reversed into round 2 is FORWARD overall,
    which is why the two signs multiply rather than either one winning.

    contig_sub_lo is the first contig base actually used: round 2 may take a SUBRANGE of a
    round-1 scaffold, so only part of a contig can be present, and a pair position outside it
    has to be dropped rather than mapped.
    """
    out = {}
    outer = c2.get(want) if c2 else None
    if outer is None:
        # no round 2: the round-1 object IS final, so its components map directly
        for s in c1.get(want, []):
            out[s["cid"]] = (s["obeg"], s["oend"], s["cbeg"], s["orient"])
        return out

    for c in outer:
        subs = [s for s in c1.get(c["cid"], [])
                if s["oend"] >= c["cbeg"] and s["obeg"] <= c["cend"]]
        if not subs:
            continue
        for s in subs:
            # the part of this round-1 component that round 2 actually used
            lo1 = max(s["obeg"], c["cbeg"])
            hi1 = min(s["oend"], c["cend"])
            # and the corresponding contig sub-interval
            if s["orient"] == "-":
                csub_lo = s["cend"] - (hi1 - s["obeg"])
            else:
                csub_lo = s["cbeg"] + (lo1 - s["obeg"])
            off_lo, off_hi = lo1 - c["cbeg"], hi1 - c["cbeg"]
            if c["orient"] == "-":
                f_lo, f_hi = c["oend"] - off_hi, c["oend"] - off_lo
            else:
                f_lo, f_hi = c["obeg"] + off_lo, c["obeg"] + off_hi
            # orientations COMPOSE: two reversals cancel
            orient = "+" if (s["orient"] == "-") == (c["orient"] == "-") else "-"
            out[s["cid"]] = (f_lo, f_hi, csub_lo, orient)
    return out


def translate(pos, pl):
    """A contig position -> a final scaffold position, or None if outside the used subrange.

    pl is (final_lo, final_hi, contig_sub_lo, orient). The offset is measured from the first
    contig base actually used, and a '-' placement reflects it from the far end.
    """
    f_lo, f_hi, c_lo, orient = pl
    off = pos - c_lo
    if off < 0 or off > (f_hi - f_lo):
        return None
    res = (f_hi - off) if orient == "-" else (f_lo + off)
    # THE RESULT MUST LAND INSIDE THIS CONTIG'S OWN FINAL INTERVAL.
    #
    # The distance invariant below is weaker than it looks: it only proves the transform is
    # RIGID, not that it is the CORRECT rigid transform. Reflecting from the wrong anchor
    # (f_lo - off instead of f_hi - off) preserves every distance and was NOT caught -- a
    # deliberately broken variant reported 0/3 mismatches. A range check catches it, because
    # a wrong anchor puts positions outside the interval the contig actually occupies.
    if res < f_lo or res > f_hi:
        return ("OUT_OF_INTERVAL", res, f_lo, f_hi)
    return res


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pairs", required=True, help="contig-space pairs.gz from FILTER_HIC_BAM")
    p.add_argument("--round1", required=True)
    p.add_argument("--round2", default="")
    p.add_argument("--scaffolds", required=True,
                   help="comma-separated final scaffold names to extract, filtered in ONE "
                        "pass over the input")
    p.add_argument("--outdir", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--pair-type", default="UU",
                   help="keep only this pair_type (default UU; the input is already filtered "
                        "to it, so this is a guard rather than a filter)")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    c1 = parse_agp(a.round1)
    c2 = parse_agp(a.round2) if a.round2 and os.path.isfile(a.round2) else None

    scafs = [s for s in a.scaffolds.split(",") if s.strip()]
    if not scafs:
        sys.exit("ERROR: --scaffolds is empty")

    # contig -> (scaffold, placement). A contig belongs to at most one final scaffold, so a
    # single flat map serves every candidate and the stream does one lookup per mate.
    where, per_scaf = {}, {}
    for sc in scafs:
        pl = contig_placements(c1, c2, sc)
        if not pl:
            sys.stderr.write("[hic_pairs] %s: no contigs resolve to %s -- skipping\n"
                             % (a.label, sc))
            continue
        per_scaf[sc] = pl
        for cid, v in pl.items():
            if cid in where:
                sys.stderr.write("[hic_pairs] WARNING contig %s maps to both %s and %s; "
                                 "keeping the first\n" % (cid, where[cid][0], sc))
                continue
            where[cid] = (sc, v)
    if not per_scaf:
        sys.exit("ERROR: none of the requested scaffolds could be resolved from the AGP(s)")
    for sc, pl in per_scaf.items():
        span = max(v[1] for v in pl.values())
        sys.stderr.write("[hic_pairs] %s %s: %d contigs, span %d bp\n"
                         % (a.label, sc, len(pl), span))

    handles = {sc: open(os.path.join(a.outdir, "%s.%s.pairs" % (a.label, sc)), "w")
               for sc in per_scaf}
    n_in = n_out = n_skip_cross = n_skip_range = n_swap = 0
    n_dist_checked = n_dist_bad = 0

    with opener(a.pairs) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            n_in += 1
            if a.pair_type and f[7] != a.pair_type:
                continue
            w1, w2 = where.get(f[1]), where.get(f[2 + 1])
            if w1 is None or w2 is None:
                continue
            if w1[0] != w2[0]:
                # both mates must be on the SAME candidate scaffold: a pair spanning two
                # different scaffolds says nothing about a junction inside either
                n_skip_cross += 1
                continue
            sc = w1[0]
            try:
                p1, p2 = int(f[2]), int(f[4])
            except ValueError:
                continue
            t1, t2 = translate(p1, w1[1]), translate(p2, w2[1])
            if t1 is None or t2 is None:
                n_skip_range += 1
                continue
            for tt, cid in ((t1, f[1]), (t2, f[3])):
                if isinstance(tt, tuple):
                    sys.exit("ERROR: contig %s position translated to %d, outside its own "
                             "final interval [%d, %d]. The coordinate mapping is wrong -- "
                             "most likely the reflection anchor or the orientation "
                             "composition. Refusing to emit a contact map that would look "
                             "plausible and be incorrect."
                             % (cid, tt[1], tt[2], tt[3]))

            # SELF-CHECK: translation is a rigid shift or reflection, so an intra-contig
            # pair's distance is invariant. Any mismatch means the arithmetic or the
            # orientation handling is wrong.
            if f[1] == f[3]:
                n_dist_checked += 1
                if abs(t2 - t1) != abs(p2 - p1):
                    n_dist_bad += 1

            s1, s2 = f[5], f[6]
            # RE-CANONICALISE to the upper triangle in FINAL coordinates. The input is upper
            # triangle in CONTIG order, which translation does not preserve.
            if (t1, t2) > (t2, t1):
                t1, t2, s1, s2 = t2, t1, s2, s1
                n_swap += 1
            handles[sc].write("%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\n"
                              % (f[0], sc, t1, sc, t2, s1, s2, f[7]))
            n_out += 1

    for h in handles.values():
        h.close()

    # chromsizes, one per scaffold: cooler cload needs it and generating it removes any
    # dependence on the input header's contig-space #chromsize lines
    for sc, pl in per_scaf.items():
        with open(os.path.join(a.outdir, "%s.%s.chromsizes" % (a.label, sc)), "w") as out:
            out.write("%s\t%d\n" % (sc, max(v[1] for v in pl.values())))

    with open(os.path.join(a.outdir, "%s.hic_pairs_audit.tsv" % a.label), "w") as out:
        out.write("# Contig-space pairs translated into candidate scaffold coordinates.\n")
        out.write("# distance_mismatches MUST be 0, but it is a WEAK check: it proves the\n")
        out.write("#   transform is rigid, not that it is the right one. Reflecting from the\n")
        out.write("#   wrong anchor preserves every distance and passes it. The range check\n")
        out.write("#   in translate() -- every position must land inside its contig's own\n")
        out.write("#   final interval -- is what catches a wrong anchor, and it aborts.\n")
        out.write("# NEITHER catches a wrong ORIENTATION: composing the two AGP levels wrongly\n")
        out.write("#   still maps inside the right interval, just backwards within it. That is\n")
        out.write("#   verified by hand and by agp_joins.py's neighbour-consistency check, not\n")
        out.write("#   here.\n")
        out.write("# swapped counts pairs re-canonicalised to the upper triangle in FINAL\n")
        out.write("#   coordinates. The input is upper triangle in CONTIG order, which\n")
        out.write("#   translation does not preserve -- without the swap the matrix is\n")
        out.write("#   silently asymmetric.\n")
        out.write("metric\tvalue\n")
        out.write("label\t%s\n" % a.label)
        out.write("pairs_read\t%d\n" % n_in)
        out.write("pairs_written\t%d\n" % n_out)
        out.write("skipped_cross_scaffold\t%d\n" % n_skip_cross)
        out.write("skipped_outside_used_subrange\t%d\n" % n_skip_range)
        out.write("swapped\t%d\n" % n_swap)
        out.write("distance_checked\t%d\n" % n_dist_checked)
        out.write("distance_mismatches\t%d\n" % n_dist_bad)
        for sc, pl in per_scaf.items():
            out.write("contigs.%s\t%d\n" % (sc, len(pl)))
            out.write("span.%s\t%d\n" % (sc, max(v[1] for v in pl.values())))

    sys.stderr.write("[hic_pairs] %s: read %d, wrote %d, swapped %d, "
                     "distance mismatches %d/%d\n"
                     % (a.label, n_in, n_out, n_swap, n_dist_bad, n_dist_checked))
    if n_dist_bad:
        sys.exit("ERROR: %d of %d intra-contig pairs changed distance under translation. The "
                 "coordinate mapping is wrong -- refusing to emit a contact map that would "
                 "look plausible and be incorrect." % (n_dist_bad, n_dist_checked))
    if not n_out:
        sys.exit("ERROR: no pairs survived. Either no contig of the requested scaffolds "
                 "appears in the pairs file, or the contig names do not match the AGP "
                 "(pairs use e.g. h1tg000001l_1).")


if __name__ == "__main__":
    main()
