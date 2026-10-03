#!/usr/bin/env python3
"""Candidate-scoped, alignment-derived transitions and exact current AGP gaps."""
import argparse
import csv
import gzip
from chimera_intervals import intervals, compatible_gaps, verdict_for_interval, arm_context

FIELDS = ["assembly", "scaffold", "name", "cut_bp", "left_chrom", "right_chrom",
          "left_component", "right_component", "n_components", "n_transitions",
          "agp_join_bp", "agp_join_distance", "agp_source", "gap_len", "callable",
          "reason", "span_bp", "vote", "candidate_verdict", "transition_lo",
          "transition_hi", "evidence_only", "chromosome_member", "location_status",
          "left_anchor_start", "left_anchor_end", "right_anchor_start", "right_anchor_end",
          "minimum_anchor_bp", "compatible_gap_count", "left_arm_aligned_bp",
          "right_arm_aligned_bp", "left_arm_fraction", "right_arm_fraction",
          "assigned_union_bp", "arm_count", "arm_pattern", "structural_status",
          "decision_scope"]

def read_rows(path):
    with open(path) as handle:
        return list(csv.DictReader((l for l in handle if l.strip() and not l.startswith("#")),
                                   delimiter="\t"))

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("joins", "round1", "paf", "assembly", "scaffold-map", "ref-name-map", "out"):
        p.add_argument("--"+key, required=True)
    p.add_argument("--round2")
    p.add_argument("--min-block", type=int, default=2000)
    p.add_argument("--component-min-bp", type=int, default=100000)
    p.add_argument("--component-margin", type=float, default=2.0)
    p.add_argument("--max-join-distance", type=int, default=250000)
    a = p.parse_args()
    if a.min_block < 1 or a.component_min_bp < 1 or a.max_join_distance < 0:
        p.error("Alignment support must be positive and interval limit nonnegative")
    all_candidates = [r for r in read_rows(a.scaffold_map) if r["assembly"] == a.assembly]
    if any("chromosome_member" not in r for r in all_candidates):
        raise ValueError("Regenerate harmonization: chromosome membership is required")
    candidates = [r for r in all_candidates if r["chromosome_member"] == "yes"
                  and r["verdict"] in ("BREAK_CANDIDATE", "REVIEW")]
    refmap = {}
    for r in read_rows(a.ref_name_map):
        n = r["new_name"]
        if n.startswith("chr") and "+" not in n:
            refmap[r["old_name"]] = n.rsplit("_", 1)[0]
    wanted = {r["scaffold"] for r in candidates}
    alns = {s: [] for s in wanted}
    opener = gzip.open if a.paf.endswith(".gz") else open
    with opener(a.paf, "rt") as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if len(f) < 12:
                raise ValueError("Malformed reference PAF")
            if f[0] not in wanted or int(f[9]) < a.min_block:
                continue
            ch = refmap.get(f[5].split("#")[-1])
            if ch:
                alns[f[0]].append((int(f[2]), int(f[3]), ch))
    gaps = read_rows(a.joins)
    output = []
    for candidate in candidates:
        sc = candidate["scaffold"]
        # Relative cap preserves support on small inferred chromosomes.
        minimum = max(1, min(a.component_min_bp, int(int(candidate["span_bp"])*0.01)))
        changes = intervals(alns[sc], minimum)
        available = [(int(g["final_cut"]), g) for g in gaps
                     if g.get("assembly") == a.assembly and g["final_object"] == sc]
        if not changes:
            row = {k: "." for k in FIELDS}
            row.update(assembly=a.assembly, scaffold=sc, name=candidate["name"],
                       span_bp=candidate["span_bp"], callable="no", evidence_only="no",
                       chromosome_member="yes", candidate_verdict="REVIEW",
                       reason="no_supported_alignment_transition", location_status="unresolved")
            output.append(row)
        for interval in changes:
            matches = compatible_gaps(interval, available)
            wide = interval["hi"]-interval["lo"] > 2*a.max_join_distance
            conflict = len(matches) == 1 and any(lo < matches[0][0] < hi for lo, hi, ch in alns[sc])
            unique = len(matches) == 1 and not wide and not conflict
            pos = matches[0][0] if unique else (interval["lo"]+interval["hi"])//2
            gap = matches[0][1] if unique else {}
            decision = verdict_for_interval(candidate, interval, len(changes)) if unique else "REVIEW"
            row = {k: "." for k in FIELDS}
            row.update(arm_context(alns[sc], interval, minimum))
            row.update(assembly=a.assembly, scaffold=sc, name=candidate["name"], cut_bp=str(pos),
                       left_chrom=interval["left"], right_chrom=interval["right"],
                       n_transitions=str(len(changes)), agp_join_bp=str(pos) if unique else ".",
                       agp_join_distance="0" if unique else ".",
                       agp_source=gap.get("lift", gap.get("source", ".")),
                       gap_len=gap.get("gap_len", "."), callable="yes" if unique else "no",
                       reason=("unique_gap_between_unambiguous_alignment_anchors" if unique else
                               "no_gap_between_alignment_anchors" if not matches else
                               "multiple_gaps_between_alignment_anchors" if len(matches) > 1 else
                               "alignment_spans_proposed_gap" if conflict else "transition_interval_too_wide"),
                       span_bp=candidate["span_bp"], vote=candidate["vote"],
                       candidate_verdict=decision, transition_lo=str(interval["lo"]),
                       transition_hi=str(interval["hi"]), evidence_only="no" if unique else "yes",
                       chromosome_member="yes",
                       location_status="unique_supported_gap" if unique else "unresolved",
                       left_anchor_start=str(interval["left_start"]), left_anchor_end=str(interval["lo"]),
                       right_anchor_start=str(interval["hi"]), right_anchor_end=str(interval["right_end"]),
                       minimum_anchor_bp=str(minimum), compatible_gap_count=str(len(matches)),
                       structural_status="alignment_transition_requires_adjudication",
                       decision_scope="single_transition_vote" if len(changes) == 1 else
                                      "junction_vote_not_available")
            output.append(row)
    with open(a.out, "w") as handle:
        w = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(output)

if __name__ == "__main__":
    main()
