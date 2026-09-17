#!/usr/bin/env python3
# ======================================================================================
# agp_joins.py
#
# Every scaffolding join, with an exact position in the FINAL scaffold's coordinates.
#
# WHY THE AGP RATHER THAN INFERENCE
# ---------------------------------
# The first design estimated a chimeric junction by binning the reference PAF at 1 Mb and
# finding where the dominant chromosome changed. Measured against the AGP on the two real
# candidates, that was wrong by 4.5 Mb and 2.0 Mb -- in both cases it landed on a
# subtelomeric region where alignment degrades and Hi-C contact is depleted, not on the glue.
#
# The AGP does not estimate anything. It records the join YaHS made, with a 100 bp
# `proximity_ligation` gap at that position by construction:
#
#   Sde-CTlk_104_hap1 scaffold_1 = scaffold_12 (41.9 Mb) + GAP + scaffold_2 (69.7 Mb)
#   Sde-CTlk_104_hap2 scaffold_3 = scaffold_36 (3.0 Mb) + GAP + scaffold_6 (43.5 Mb, REVERSE,
#                                  subrange 4.3-47.8 Mb) + GAP + scaffold_16 (30.2 Mb)
#
# and the PAF then says which reference chromosome each COMPONENT belongs to, which is a much
# better-posed question than asking where a boundary falls.
#
# THESE ARE CANDIDATES, NOT THINGS TO BREAK
# -----------------------------------------
# Sde-CPla_115_hap1 has 722 round-1 joins and 197 round-2 joins. Breaking at every AGP join
# would shatter the assembly back to contigs and undo scaffolding entirely. The AGP supplies a
# SUPERSET of candidate positions; the concordance vote remains the gate, and the span and
# footprint floors exclude fragment-scale joins.
#
# CHAINING IS REQUIRED, NOT OPTIONAL
# ----------------------------------
# Round 1 makes 1,830-1,891 joins per haplotype here; round 2 makes 32-54. So the great
# majority of joins -- and of potential chimeras -- are round-1, recorded only in the round-1
# AGP in round-1 coordinates. Our two known candidates are round-2, which is luck rather than
# the general case. Without chaining, ~97% of the search space is invisible.
#
# Chaining a round-1 join into round-2 coordinates has to handle three things the naive
# version gets wrong: components entering as a SUBRANGE (scaffold_6 enters as 4306001-47785970,
# so YaHS broke it and used the middle), REVERSE orientation (scaffold_6 again), and round-1
# scaffolds that never entered round 2 at all.
#
# Dependencies: python3 stdlib.
#
# USAGE
#   agp_joins.py --round1 r1.agp [--round2 r2.agp] --assembly <id> --out joins.tsv
#   agp_joins.py --round1 r1.agp --round2 r2.agp --assembly <id> --out joins.tsv \
#       --require-gap-type proximity_ligation --require-gap-len 100
# ======================================================================================

import argparse
import os
import sys


# --------------------------------------------------------------------------------------
def parse_agp(path):
    """AGP -> (components, gaps).

    components[object] = [ {part, obeg, oend, cid, cbeg, cend, orient}, ... ]
    gaps[object]       = [ {part, obeg, oend, length, gap_type, evidence}, ... ]

    AGP columns: object object_beg object_end part_number component_type
                 then either  component_id component_beg component_end orientation   (W/D/etc)
                 or           gap_length gap_type linkage linkage_evidence            (N/U)

    All coordinates are 1-based inclusive, as AGP specifies. Converted to 0-based only at
    the point of emitting a cut position, and that conversion is stated where it happens --
    an off-by-one here silently moves every cut.
    """
    comps, gaps = {}, {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            obj, ob, oe, part, ctype = f[0], int(f[1]), int(f[2]), int(f[3]), f[4]
            if ctype in ("N", "U"):
                gaps.setdefault(obj, []).append({
                    "part": part, "obeg": ob, "oend": oe,
                    "length": int(f[5]) if f[5].isdigit() else (oe - ob + 1),
                    "gap_type": f[6], "evidence": f[8] if len(f) > 8 else "."})
            else:
                comps.setdefault(obj, []).append({
                    "part": part, "obeg": ob, "oend": oe, "cid": f[5],
                    "cbeg": int(f[6]), "cend": int(f[7]), "orient": f[8]})
    for d in (comps, gaps):
        for k in d:
            d[k].sort(key=lambda x: x["part"])
    return comps, gaps


# --------------------------------------------------------------------------------------
def build_lift(comps2):
    """round-1 component id -> the round-2 placements that used it.

    A round-1 scaffold can appear MORE THAN ONCE if round 2 broke it and placed the pieces
    separately, so this maps to a list rather than a single placement. Each entry carries the
    subrange actually used, which is why a round-1 position has to be checked for membership
    rather than assumed present.
    """
    lift = {}
    for obj, parts in comps2.items():
        for c in parts:
            lift.setdefault(c["cid"], []).append({
                "obj": obj, "obeg": c["obeg"], "oend": c["oend"],
                "cbeg": c["cbeg"], "cend": c["cend"], "orient": c["orient"]})
    return lift


def lift_position(lift, r1_scaffold, r1_pos):
    """A position on a round-1 scaffold -> (round-2 object, position), or None.

    r1_pos is 1-based, as AGP coordinates are.

    FORWARD: the component occupies object [obeg, oend] and corresponds to round-1
             [cbeg, cend], so offset within the component is preserved:
                 obj_pos = obeg + (r1_pos - cbeg)

    REVERSE: the component is reverse-complemented into the object, so the offset is
             measured from the other end:
                 obj_pos = oend - (r1_pos - cbeg)
             This is the case that a careful reading gets wrong and a test catches --
             scaffold_6 enters Sde-CTlk_104_hap2 scaffold_3 reversed, in a real candidate.

    Returns None when the position falls outside every subrange that was used: round 2 may
    have broken the round-1 scaffold and discarded the piece containing it.
    """
    out = []
    for p in lift.get(r1_scaffold, []):
        if not (p["cbeg"] <= r1_pos <= p["cend"]):
            continue
        off = r1_pos - p["cbeg"]
        if p["orient"] == "-":
            obj_pos = p["oend"] - off
        else:
            obj_pos = p["obeg"] + off
        out.append((p["obj"], obj_pos, p["orient"]))
    return out


# --------------------------------------------------------------------------------------
def joins_from_agp(comps, gaps, source):
    """Every gap that JOINS two components -> a join record, in this AGP's own coordinates.

    A gap only counts as a join when it has a component on both sides. A leading or trailing
    gap is padding, not a join, and breaking there would produce an empty piece.

    The cut position is the gap's MIDPOINT: cutting inside the gap severs no real sequence,
    which is the whole reason the AGP position is preferable to an inferred one.
    """
    out = []
    for obj, gl in gaps.items():
        cs = comps.get(obj, [])
        if len(cs) < 2:
            continue
        by_part = {c["part"]: c for c in cs}
        for g in gl:
            left = by_part.get(g["part"] - 1)
            right = by_part.get(g["part"] + 1)
            if not (left and right):
                continue
            out.append({
                "source": source,
                "object": obj,
                "gap_part": g["part"],
                "gap_beg": g["obeg"],
                "gap_end": g["oend"],
                "gap_len": g["length"],
                "gap_type": g["gap_type"],
                "evidence": g["evidence"],
                "cut_pos": (g["obeg"] + g["oend"]) // 2,
                "left_cid": left["cid"], "left_bp": left["oend"] - left["obeg"] + 1,
                "left_orient": left["orient"],
                "right_cid": right["cid"], "right_bp": right["oend"] - right["obeg"] + 1,
                "right_orient": right["orient"],
            })
    return out


# --------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--round1", required=True)
    p.add_argument("--round2", default="",
                   help="omit when round 2 did not run; round-1 objects are then final")
    p.add_argument("--assembly", required=True)
    p.add_argument("--out", required=True)
    # AGP gap columns are: 6=gap_length 7=gap_type 8=linkage 9=linkage_evidence. A YaHS
    # scaffolding join is gap_type `scaffold` with evidence `proximity_ligation` -- those are
    # DIFFERENT fields, and conflating them made the filter reject every join.
    p.add_argument("--require-gap-type", default="",
                   help="keep only joins with this gap_type (column 7; YaHS uses 'scaffold')")
    p.add_argument("--require-evidence", default="",
                   help="keep only joins with this linkage_evidence (column 9). "
                        "'proximity_ligation' selects Hi-C scaffolding joins. Measured here: "
                        "all 32 round-2 gaps are `scaffold proximity_ligation`, 100 bp -- read "
                        "from the AGP rather than assumed, so a different scaffolder simply "
                        "yields no matches and every candidate goes to REVIEW.")
    p.add_argument("--require-gap-len", type=int, default=0,
                   help="keep only joins with exactly this gap length (0 = any)")
    a = p.parse_args()

    if not os.path.isfile(a.round1):
        sys.exit("ERROR: no such file: %s" % a.round1)
    c1, g1 = parse_agp(a.round1)
    j1 = joins_from_agp(c1, g1, "round1")
    sys.stderr.write("[agp_joins] %s: round1 %d objects, %d joins\n"
                     % (a.assembly, len(c1), len(j1)))

    rows = []
    if a.round2 and os.path.isfile(a.round2):
        c2, g2 = parse_agp(a.round2)
        j2 = joins_from_agp(c2, g2, "round2")
        lift = build_lift(c2)
        sys.stderr.write("[agp_joins] %s: round2 %d objects, %d joins\n"
                         % (a.assembly, len(c2), len(j2)))

        # round-2 joins are already in final-object coordinates
        for j in j2:
            j["final_object"] = j["object"]
            j["final_cut"] = j["cut_pos"]
            j["lift"] = "native"
            rows.append(j)

        # round-1 joins have to be carried through. Round 1 dominates the join count
        # (1,830-1,891 vs 32-54 per haplotype here), so skipping this would hide ~97% of
        # the candidate space.
        n_lift, n_drop, n_multi = 0, 0, 0
        for j in j1:
            hits = lift_position(lift, j["object"], j["cut_pos"])
            if not hits:
                n_drop += 1
                continue
            if len(hits) > 1:
                n_multi += 1
            for obj, pos, orient in hits:
                k = dict(j)
                k["final_object"] = obj
                k["final_cut"] = pos
                k["lift"] = "round1_via_round2%s" % ("_rev" if orient == "-" else "")
                if orient == "-":
                    # SWAP. left/right are recorded in ROUND-1 order, but a reverse-oriented
                    # component is reverse-complemented into the object, so the round-1 left
                    # neighbour ends up at the HIGHER final coordinate. Visible in the real
                    # output: consecutive _rev joins listed h2tg001770l as the left of one and
                    # the right of the next, which cannot both hold.
                    #
                    # This is load-bearing rather than cosmetic: the next step asks which
                    # reference chromosome lies left versus right of each join, and unswapped
                    # labels invert that answer for every join inside a reverse component --
                    # 100+ of them on the Sde-CTlk_104_hap2 candidate alone.
                    for lo, hi in (("left_cid", "right_cid"), ("left_bp", "right_bp"),
                                   ("left_orient", "right_orient")):
                        k[lo], k[hi] = k[hi], k[lo]
                rows.append(k)
                n_lift += 1
        sys.stderr.write("[agp_joins] %s: lifted %d round1 joins (%d dropped -- the round-1 "
                         "scaffold or that part of it did not enter round 2; %d placed more "
                         "than once)\n" % (a.assembly, n_lift, n_drop, n_multi))
    else:
        for j in j1:
            j["final_object"] = j["object"]
            j["final_cut"] = j["cut_pos"]
            j["lift"] = "round1_final"
            rows.append(j)
        sys.stderr.write("[agp_joins] %s: no round2 AGP; round1 objects are final\n"
                         % a.assembly)

    kept = []
    for r in rows:
        if a.require_gap_type and r["gap_type"] != a.require_gap_type:
            continue
        if a.require_evidence and r["evidence"] != a.require_evidence:
            continue
        if a.require_gap_len and r["gap_len"] != a.require_gap_len:
            continue
        kept.append(r)
    if len(kept) != len(rows):
        sys.stderr.write("[agp_joins] %s: %d of %d joins pass the gap filters\n"
                         % (a.assembly, len(kept), len(rows)))

    cols = ["assembly", "final_object", "final_cut", "source", "lift", "object", "cut_pos",
            "gap_beg", "gap_end", "gap_len", "gap_type", "evidence",
            "left_cid", "left_bp", "left_orient", "right_cid", "right_bp", "right_orient"]
    with open(a.out, "w") as out:
        out.write("# Every scaffolding join, positioned in the FINAL scaffold's coordinates.\n")
        out.write("# CANDIDATES ONLY. Breaking at every join would undo scaffolding -- one\n")
        out.write("#   haplotype here has 919 of them. The concordance vote is the gate; this\n")
        out.write("#   file only says WHERE a join is, exactly, rather than estimating it.\n")
        out.write("# final_cut is 1-based, the gap MIDPOINT: cutting inside the gap severs no\n")
        out.write("#   real sequence. It is a round-2-coordinate position -- gap filling and\n")
        out.write("#   teloclip shift it slightly afterwards, so the consumer snaps to the\n")
        out.write("#   nearest gap of gap_len in the current assembly. Measured offsets on\n")
        out.write("#   the known candidates: -9 kb, +0 kb, -0 kb.\n")
        out.write("# lift: native = a round-2 join. round1_via_round2 = carried through, and\n")
        out.write("#   _rev means the round-1 scaffold entered reversed: the offset was\n")
        out.write("#   measured from the far end AND left_cid/right_cid were swapped, because\n")
        out.write("#   reverse-complementing puts the round-1 left neighbour at the higher\n")
        out.write("#   final coordinate. left_cid is always the lower-coordinate side.\n")
        out.write("\t".join(cols) + "\n")
        for r in sorted(kept, key=lambda x: (x["final_object"], x["final_cut"])):
            r["assembly"] = a.assembly
            out.write("\t".join(str(r.get(c, ".")) for c in cols) + "\n")
    sys.stderr.write("[agp_joins] %s: wrote %d joins -> %s\n"
                     % (a.assembly, len(kept), os.path.basename(a.out)))


if __name__ == "__main__":
    main()
