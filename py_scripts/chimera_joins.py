#!/usr/bin/env python3
# ======================================================================================
# chimera_joins.py
#
# Which scaffolding joins separate two different reference chromosomes -- the chimeric ones.
#
# HOW THIS WORKS, AND WHY EACH PIECE IS NEEDED
# --------------------------------------------
#   1. AGP chain          exact join positions, and the component intervals between them
#   2. PAF per component  which reference chromosome each component belongs to
#   3. fill               a component with no alignment is not evidence of a change
#   4. transitions        adjacent components with DIFFERENT chromosomes -> a chimeric join
#
# All four are load-bearing. Measured on the two real candidates:
#
#   Sde-CTlk_104_hap1 chr5_1+chr9_1, 213 components, 1 real transition at 37,144,373
#   Sde-CTlk_104_hap2 chr6_3+chr12_1, 174 components, 1 real transition at 43,694,728
#
# Both are ROUND-1 joins, so the chaining in agp_joins.py is required: round 1 makes
# 1,830-1,891 joins per haplotype against round 2's 32-54, and without lifting them ~97% of
# the search space -- including both real junctions -- is invisible.
#
# WHY NOT THE AGP ALONE, AND WHY NOT INFERENCE ALONE
# --------------------------------------------------
# The AGP says where joins ARE but not which ones matter: Sde-CTlk_104_hap2 scaffold_3 has
# 173 joins and one is chimeric. Cutting at all of them would undo scaffolding.
#
# Inference alone lands near but not on the junction. Estimates for chr5_1+chr9_1:
#
#   PAF binning at 1 Mb      37.0 Mb    -144 kb from the truth
#   Hi-C cross-contact min   37.4 Mb    +256 kb
#   the actual round-1 join  37,144,322
#
# Close, but not a position you can cut at -- and both sit in sequence, not in the 100 bp gap
# that the join actually is. The combination is what works: inference cannot give an exact
# coordinate and the AGP cannot say which join matters.
#
# (I spent a round arguing the junction was the ROUND-2 join at 41,919,050, 4.77 Mb away, and
# built an explanation about repeat-rich regions misleading the PAF to account for a
# discrepancy I had manufactured. Every component from 37.14 to 41.92 Mb is chr9, so that
# join is chr9->chr9 and not chimeric. Hence this script tests every join rather than
# assuming the largest or the most recent is the relevant one.)
#
# THREE RULES THAT LOOK LIKE DETAILS AND ARE NOT
# ----------------------------------------------
#   `unplaced` is NOT a chromosome. It is the reference's own unassigned sequence, so it says
#   nothing about which chromosome a component belongs to. Treating it as one produced a
#   spurious out-and-back transition pair on an 86 kb component. Only `^chr` targets vote;
#   measured, `unplaced` is the only non-chr target and totals 0.2 Mb per assembly.
#
#   Unassigned components are FILLED from their neighbours. h2tg000137l_1 aligns over 0.83 Mb
#   of its 4.85 Mb span, and components near a junction have less -- unfilled, every
#   repeat-rich contig reads as a transition.
#
#   A component needs a minimum aligned bp AND a dominance margin before it votes. Without
#   the margin, a component with 60 kb on one chromosome and 55 kb on another would pick a
#   side arbitrarily.
#
# Dependencies: python3 stdlib.
#
# USAGE
#   chimera_joins.py --joins <agp_joins.tsv> --round1 r1.agp [--round2 r2.agp] \
#       --paf ref_vs_asm.paf.gz --assembly <id> --scaffold-map <candidates.tsv> \
#       --out joins_called.tsv
# ======================================================================================

import argparse
import gzip
import os
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_rows(path, sep="\t"):
    rows, hdr = [], None
    if not path or not os.path.isfile(path):
        return rows
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split(sep)
            if hdr is None:
                hdr = f
                continue
            if len(f) == len(hdr):
                rows.append(dict(zip(hdr, f)))
    return rows


# --------------------------------------------------------------------------------------
def parse_agp_components(path):
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


def final_components(c1, c2, obj):
    """Component intervals of one FINAL scaffold, at contig granularity.

    Each round-2 component is a round-1 scaffold or a subrange of one; it is expanded into
    that scaffold's own components so the unit is a contig. A contig cannot straddle a
    scaffolding join by construction, which is what makes it the right unit for chromosome
    assignment -- a fixed window around a join can and does straddle one.

    The reverse case is the one to get right: a component placed with orient '-' is
    reverse-complemented into the object, so its round-1 sub-intervals map to final
    coordinates in REVERSE order and measured from the far end.
    """
    out = []
    for c in c2.get(obj, []) if c2 else [dict(part=1, obeg=None, oend=None, cid=obj,
                                              cbeg=None, cend=None, orient="+")]:
        subs = c1.get(c["cid"], [])
        if c["obeg"] is None:                       # no round 2: round-1 object is final
            for s in subs:
                out.append((s["obeg"], s["oend"], s["cid"], "+"))
            continue
        use = [s for s in subs if s["oend"] >= c["cbeg"] and s["obeg"] <= c["cend"]]
        if not use:
            out.append((c["obeg"], c["oend"], c["cid"], c["orient"]))
            continue
        for s in use:
            lo = max(s["obeg"], c["cbeg"])
            hi = min(s["oend"], c["cend"])
            o_lo, o_hi = lo - c["cbeg"], hi - c["cbeg"]
            if c["orient"] == "-":
                f_lo, f_hi = c["oend"] - o_hi, c["oend"] - o_lo
            else:
                f_lo, f_hi = c["obeg"] + o_lo, c["obeg"] + o_hi
            out.append((f_lo, f_hi, s["cid"], c["orient"]))
    out.sort()
    return out


# --------------------------------------------------------------------------------------
def read_ref_name_map(path):
    """The REFERENCE's own harmonized_name_map.tsv -> {reference scaffold: consensus chrN}.

    THIS TRANSLATION IS THE POINT, and omitting it was a real error rather than a detail.
    Harmonization maintains two namespaces: refN is a reference-frame PIECE, chrN a CONSENSUS
    chromosome built from the join graph's connected components across voters. A consensus
    chromosome can span several reference pieces when the voters agree those pieces belong
    together -- so a scaffold aligning to ref5 on one side and ref12 on the other is
    CORRECTLY JOINING TWO PIECES OF chr5, not chimeric. Assigning components to reference
    pieces would call that a chimera and destroy real assembly work.
    (Measured on this cohort every consensus chromosome is exactly one reference piece, so
    the translation is currently the identity -- but it is not on a fragmented reference,
    which is the case the consensus frame exists to handle.)

    It also matters mechanically. The PAF source determines the target names:
      harmonization's own PAFs align the INPUT fastas, so targets are `scaffold_N`
      PAIRWISE_ALIGNMENT's PAFs align FINALIZED assemblies, so targets are `chr5_1`
    Only the first is guaranteed to exist on a detection-only run, and a `^chr` test on its
    targets matches nothing -- silently, because zero transitions is a legitimate result.

    A target absent from the map, or mapping to `unplaced*`, is unassigned: unplaced sequence
    is the reference's own unassigned material and says nothing about which chromosome a
    component belongs to.
    """
    out = {}
    for r in read_rows(path):
        old, new = r.get("old_name", ""), r.get("new_name", "")
        if not old or not new or new.startswith("unplaced"):
            continue
        if "+" in new:
            # a composite IN THE REFERENCE spans two chromosomes, so it cannot name one
            continue
        out[old] = new.rsplit("_", 1)[0] if "_" in new else new
    return out


def paf_by_query(paf, wanted, min_block, ref_map):
    """query -> [(qstart, qend, consensus_chrom, block), ...].

    Targets are translated through the reference name map, so this works for either PAF
    source and yields consensus chromosomes either way. An untranslatable target is dropped
    rather than guessed at.
    """
    out = {}
    n_seen, n_kept = 0, 0
    with opener(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            q = f[0]
            if q not in wanted:
                continue
            try:
                blk = int(f[9])
                if blk < min_block:
                    continue
                qs, qe = int(f[2]), int(f[3])
            except ValueError:
                continue
            n_seen += 1
            tgt = f[5].split("#")[-1]
            ch = ref_map.get(tgt)
            if ch is None:
                # already a consensus name? PAIRWISE_ALIGNMENT PAFs look like chr5_1
                if tgt.startswith("chr"):
                    ch = tgt.rsplit("_", 1)[0] if "_" in tgt else tgt
                else:
                    continue
            n_kept += 1
            out.setdefault(q, []).append((qs, qe, ch, blk))
    if n_seen and not n_kept:
        sys.exit("ERROR: %d PAF records matched the candidate scaffolds but NONE had a "
                 "translatable target. The reference name map probably does not match this "
                 "PAF's target namespace -- harmonization's PAFs use `scaffold_N`, "
                 "PAIRWISE_ALIGNMENT's use `chr5_1`. Failing rather than reporting zero "
                 "transitions, which would look like a clean result." % n_seen)
    return out


def assign(comps, aln, min_bp, margin):
    """A reference chromosome per component, then filled from neighbours.

    Returns [ {lo, hi, cid, orient, called, filled, top_bp, second_bp}, ... ]
    """
    rows = []
    for lo, hi, cid, orient in comps:
        per = {}
        for qs, qe, ch, blk in aln:
            if qe > lo and qs < hi:                 # 0/1-based mixing is immaterial here:
                per[ch] = per.get(ch, 0) + blk      # an overlap test, not a coordinate
        top = sorted(per.items(), key=lambda x: -x[1])
        called, t0, t1 = ".", (top[0][1] if top else 0), (top[1][1] if len(top) > 1 else 0)
        if top and t0 >= min_bp and (len(top) == 1 or t0 >= margin * t1):
            called = top[0][0]
        rows.append(dict(lo=lo, hi=hi, cid=cid, orient=orient, called=called,
                         top_bp=t0, second_bp=t1, per=per))

    # ---- BLOCK assignment before filling -------------------------------------------
    # The per-component threshold is right for a repeat-rich contig in the middle of a
    # chromosome, but wrong for a RUN of small components that collectively carry plenty of
    # alignment. Sde-CTlk_104_hap2 scaffold_3 begins with scaffold_36: 3.0 Mb of contigs none
    # of which clears 100 kb alone, but 1,311 kb aligned to a DIFFERENT reference scaffold
    # than the rest. Filled from its neighbour it read as the same chromosome and its junction
    # -- an independently Hi-C-supported one, cross-contact 0.534, MORE depleted than the cut
    # that was found -- disappeared.
    #
    # So consecutive unassigned components are pooled and the POOL is tested against the same
    # threshold. That keeps the per-component guard everywhere it matters and stops a run of
    # small pieces being absorbed into whatever happens to border it.
    i = 0
    while i < len(rows):
        if rows[i]["called"] != ".":
            i += 1
            continue
        j = i
        while j < len(rows) and rows[j]["called"] == ".":
            j += 1
        blk = rows[i:j]
        per = {}
        for r in blk:
            for ch, bp in (r.get("per") or {}).items():
                per[ch] = per.get(ch, 0) + bp
        top = sorted(per.items(), key=lambda x: -x[1])
        if top and top[0][1] >= min_bp and (len(top) == 1 or top[0][1] >= margin * top[1][1]):
            for r in blk:
                r["called"] = top[0][0]
                r["block_called"] = "%s(%d bp over %d components)" % (top[0][0], top[0][1],
                                                                     len(blk))
        i = j

    known = [i for i, r in enumerate(rows) if r["called"] != "."]
    for i, r in enumerate(rows):
        if r["called"] != ".":
            r["filled"] = r["called"]
        elif known:
            r["filled"] = rows[min(known, key=lambda k: abs(k - i))]["called"]
        else:
            r["filled"] = "."
    return rows


def transitions(rows):
    out = []
    for a, b in zip(rows, rows[1:]):
        if a["filled"] != b["filled"] and "." not in (a["filled"], b["filled"]):
            out.append((b["lo"], a["filled"], b["filled"], a["cid"], b["cid"],
                        a["hi"], b["lo"]))
    return out


# --------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--joins", required=True, help="agp_joins.py output")
    p.add_argument("--round1", required=True)
    p.add_argument("--round2", default="")
    p.add_argument("--paf", required=True, help="reference_vs_assembly PAF (.gz ok)")
    p.add_argument("--assembly", required=True)
    p.add_argument("--scaffold-map", required=True,
                   help="chimera_candidates.tsv: supplies which scaffolds to test and the "
                        "harmonized name the PAF uses as its query")
    p.add_argument("--out", required=True)
    p.add_argument("--ref-name-map", required=True,
                   help="the REFERENCE's harmonized_name_map.tsv. Translates PAF targets to "
                        "CONSENSUS chromosomes, which is what a transition must be measured "
                        "in: a consensus chromosome can span several reference pieces, and a "
                        "scaffold joining two pieces of one chromosome is not chimeric.")
    p.add_argument("--min-block", type=int, default=2000,
                   help="ignore PAF records below this block length (default 2 kb)")
    p.add_argument("--component-min-bp", type=int, default=100000,
                   help="a component needs this much aligned bp to vote (default 100 kb)")
    p.add_argument("--component-margin", type=float, default=2.0,
                   help="and this dominance ratio over the runner-up (default 2.0), so a "
                        "component split near-evenly between two chromosomes abstains "
                        "rather than picking arbitrarily")
    p.add_argument("--max-join-distance", type=int, default=250000,
                   help="a transition must have an AGP join within this distance to be "
                        "callable (default 250 kb). Further away means the chromosome change "
                        "is INSIDE a contig -- a contig-level mis-assembly, which is "
                        "Inspector's domain and cannot be cut at a scaffolding gap.")
    a = p.parse_args()

    cand = [r for r in read_rows(a.scaffold_map) if r.get("assembly") == a.assembly]
    if not cand:
        sys.stderr.write("[chimera_joins] %s: no candidate scaffolds; nothing to test\n"
                         % a.assembly)
    # EVERY candidate row is tested, whatever its verdict -- detection reports, the gate
    # decides. Sde-CPla_115_hap1 yields 79 chimeric joins across ~80 scaffolds, and that is a
    # real assembly-quality finding worth surfacing rather than suppressing at detection. But
    # the verdict and span travel with each row so BREAK_CHIMERAS can act without re-deriving
    # them.
    scafs = {r["scaffold"]: r.get("name", r["scaffold"]) for r in cand}
    meta = {r["scaffold"]: r for r in cand}
    # BOTH namespaces are accepted, because the PAF source decides which appears:
    #   harmonization's PAFs align the INPUT fastas  -> query `scaffold_1`, target `scaffold_5`
    #   PAIRWISE_ALIGNMENT aligns FINALIZED assemblies -> query `chr5_1+chr9_1`, target `chr5_1`
    # The targets were already translated through the reference name map; the QUERIES were
    # not, so a harmonization PAF matched nothing and every scaffold came back "no assignable
    # component" -- correctly reported as not tested, but for the wrong reason.
    qnames = set(scafs.keys()) | set(scafs.values())
    # and back to the scaffold, whichever name the PAF used
    q_to_scaf = {}
    for s, n in scafs.items():
        q_to_scaf[s] = s
        q_to_scaf[n] = s

    c1 = parse_agp_components(a.round1)
    c2 = parse_agp_components(a.round2) if a.round2 and os.path.isfile(a.round2) else None
    ref_map = read_ref_name_map(a.ref_name_map)
    if not ref_map:
        sys.exit("ERROR: no usable rows in %s -- without the reference name map, PAF targets "
                 "cannot be resolved to consensus chromosomes." % a.ref_name_map)
    sys.stderr.write("[chimera_joins] reference name map: %d scaffolds -> %d consensus "
                     "chromosomes\n" % (len(ref_map), len(set(ref_map.values()))))
    aln = paf_by_query(a.paf, qnames, a.min_block, ref_map)
    if cand and not aln:
        sys.exit("ERROR: the PAF has no records for ANY of the %d candidate scaffolds. "
                 "Tried both namespaces: %s ... and %s ... . The PAF's query names are in a "
                 "third namespace, or it is the wrong PAF for this assembly. Failing rather "
                 "than reporting every scaffold as untested."
                 % (len(scafs), sorted(scafs.keys())[:3], sorted(scafs.values())[:3]))
    joins = [r for r in read_rows(a.joins) if r.get("assembly") == a.assembly]

    n_unassigned = 0
    cols = ["assembly", "scaffold", "name", "cut_bp", "left_chrom", "right_chrom",
            "left_component", "right_component", "n_components", "n_transitions",
            "agp_join_bp", "agp_join_distance", "agp_source", "gap_len", "callable", "reason",
            # carried from the candidates file so BREAK_CHIMERAS gates without re-deriving
            "span_bp", "vote", "candidate_verdict"]
    n_call = 0
    with open(a.out, "w") as out:
        out.write("# Joins that separate two DIFFERENT reference chromosomes.\n")
        out.write("# cut_bp is an AGP join position -- the midpoint of a 100 bp scaffolding\n")
        out.write("#   gap -- so cutting there severs no real sequence. Inference gets close\n")
        out.write("#   (-144 kb and +256 kb on the known case) but lands in sequence.\n")
        out.write("# callable=no with a large agp_join_distance means the chromosome changes\n")
        out.write("#   INSIDE a contig: a contig-level mis-assembly, not a scaffolding one.\n")
        out.write("#   That cannot be cut at a gap and belongs to Inspector.\n")
        out.write("# THIS IS NOT A DECISION. candidate_verdict, carried from\n")
        out.write("#   chimera_candidates.tsv, is the gate; this says WHICH join and WHERE.\n")
        out.write("# Only BREAK_CANDIDATE rows are cut. A REVIEW verdict means the vote could\n")
        out.write("#   not separate an artifact from real biology -- in particular a\n")
        out.write("#   POLYMORPHIC fusion, where the sister haplotype also carries the\n")
        out.write("#   junction. Those are never broken automatically.\n")
        out.write("# Every candidate is reported regardless of verdict: Sde-CPla_115_hap1 has\n")
        out.write("#   79 chimeric joins across ~80 scaffolds, which is an assembly-quality\n")
        out.write("#   finding in its own right and should not be hidden by the span floor.\n")
        out.write("\t".join(cols) + "\n")

        for scaf, qname in sorted(scafs.items()):
            comps = final_components(c1, c2, scaf)
            if not comps:
                sys.stderr.write("[chimera_joins] %s %s: no components in the AGP\n"
                                 % (a.assembly, scaf))
                continue
            # the PAF may key on either name; take whichever it used
            al = aln.get(scaf) or aln.get(qname) or []
            rows = assign(comps, al, a.component_min_bp, a.component_margin)
            if not any(r["called"] != "." for r in rows):
                # zero transitions is a legitimate result, so an unassignable scaffold must
                # not be reported as one -- it would read as "no chimeras found".
                sys.stderr.write("[chimera_joins] %s %s (%s): WARNING no component could be "
                                 "assigned a chromosome from %d alignment(s). Either the PAF "
                                 "has no records for %r or %r, or every component is below "
                                 "--component-min-bp.\n"
                                 % (a.assembly, scaf, qname, len(al), scaf, qname))
                n_unassigned += 1
                continue
            tr = transitions(rows)
            sj = sorted([j for j in joins if j.get("final_object") == scaf],
                        key=lambda x: int(x["final_cut"]))
            sys.stderr.write("[chimera_joins] %s %s (%s): %d components, %d transition(s), "
                             "%d AGP joins\n"
                             % (a.assembly, scaf, qname, len(rows), len(tr), len(sj)))

            for pos, lch, rch, lcid, rcid, lhi, rlo in tr:
                # the AGP join nearest the transition. The transition is bracketed by two
                # components, so the join between them is the one to cut at.
                best, bd = None, None
                for j in sj:
                    d = abs(int(j["final_cut"]) - pos)
                    if bd is None or d < bd:
                        best, bd = j, d
                callable_ = "yes" if (best and bd <= a.max_join_distance) else "no"
                reason = ("" if callable_ == "yes" else
                          ("nearest AGP join %d bp away -- the change is inside a contig"
                           % bd if best else "no AGP join on this scaffold"))
                if callable_ == "yes":
                    n_call += 1
                rec = dict(assembly=a.assembly, scaffold=scaf, name=qname,
                           cut_bp=(best["final_cut"] if best else pos),
                           left_chrom=lch, right_chrom=rch,
                           left_component=lcid, right_component=rcid,
                           n_components=len(rows), n_transitions=len(tr),
                           agp_join_bp=(best["final_cut"] if best else "."),
                           agp_join_distance=(bd if best is not None else "."),
                           agp_source=(best["lift"] if best else "."),
                           gap_len=(best["gap_len"] if best else "."),
                           callable=callable_, reason=(reason or "."),
                           span_bp=meta[scaf].get("span_bp", "."),
                           vote=meta[scaf].get("vote", "."),
                           candidate_verdict=meta[scaf].get("verdict", "."))
                out.write("\t".join(str(rec[c]) for c in cols) + "\n")

    by_verdict = {}
    for r in read_rows(a.out):
        if r.get("callable") == "yes":
            k = r.get("candidate_verdict", "?")
            by_verdict[k] = by_verdict.get(k, 0) + 1
    if n_unassigned:
        sys.stderr.write("[chimera_joins] %s: %d scaffold(s) had NO assignable component -- "
                         "those were not tested, not cleared\n" % (a.assembly, n_unassigned))
    sys.stderr.write("[chimera_joins] %s: %d callable chimeric join(s) -> %s\n"
                     % (a.assembly, n_call, os.path.basename(a.out)))
    if by_verdict:
        sys.stderr.write("[chimera_joins] %s: by candidate verdict -- %s\n"
                         % (a.assembly, ", ".join("%s=%d" % kv
                                                  for kv in sorted(by_verdict.items()))))


if __name__ == "__main__":
    main()
