#!/usr/bin/env python3
"""Conservative older-join recovery using exact, unique minimap2 flank alignments.
Coordinates are 0-based half-open internally. Every input join gets an audit row.
Recovered joins are REVIEW only, never automatic cutting authorization.
"""
from chimera_schema import read_table, write_called
import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


def records(path):
    name, parts = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif name is not None:
                parts.append(line.strip())
        if name is not None:
            yield name, "".join(parts).upper()


def rows(path):
    with open(path) as handle:
        reader = csv.DictReader((s for s in handle if s.strip() and not s.startswith("#")),
                                delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def checksum(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def prepare(agp, fasta, flank, minimum, prefix):
    components, gaps = defaultdict(dict), []
    with open(agp) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) != 9:
                raise ValueError("Malformed AGP")
            obj, lo, hi, part = f[0], int(f[1])-1, int(f[2]), int(f[3])
            if f[4] in ("N", "U"):
                gaps.append(dict(source_object=obj, source_lo=lo, source_hi=hi,
                                 part=part, evidence=f[8]))
            else:
                components[obj][part] = (lo, hi)
    by_obj = defaultdict(list)
    joins = []
    for gap in gaps:
        cs = components[gap["source_object"]]
        if gap["part"]-1 not in cs or gap["part"]+1 not in cs:
            continue
        gap["join_id"] = "J%08d" % len(joins)
        gap["status"] = "awaiting_alignment"
        gap["left_lo"] = max(cs[gap["part"]-1][0], gap["source_lo"]-flank)
        gap["right_hi"] = min(cs[gap["part"]+1][1], gap["source_hi"]+flank)
        joins.append(gap)
        by_obj[gap["source_object"]].append(gap)
    seen = set()
    with open(prefix+".flanks.fa", "w") as out:
        for name, seq in records(fasta):
            if name in seen:
                raise ValueError("Duplicate FASTA name")
            seen.add(name)
            for j in by_obj.get(name, []):
                if j["right_hi"] > len(seq) or set(seq[j["source_lo"]:j["source_hi"]]) != {"N"}:
                    raise ValueError("Original AGP gap disagrees with first-round FASTA")
                left = seq[j["left_lo"]:j["source_lo"]]
                right = seq[j["source_hi"]:j["right_hi"]]
                if min(len(left), len(right)) < minimum or re.search("[^ACGT]", left+right):
                    j["status"] = "unresolved_short_or_ambiguous_source_flank"
                    continue
                for side, flankseq in (("L", left), ("R", right)):
                    out.write(">"+j["join_id"]+"_"+side+"\n"+flankseq+"\n")
    if set(by_obj)-seen:
        raise ValueError("AGP object absent from first-round FASTA")
    Path(prefix+".source.json").write_text(json.dumps(
        dict(joins=joins, source_sha256=checksum(fasta)), indent=2)+"\n")


def exact_hits(path):
    """Count all emitted full-length exact placements, including low-MAPQ secondaries."""
    hits = defaultdict(dict)
    with open(path) as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if len(f) < 12:
                raise ValueError("Malformed flank PAF")
            qlen, qs, qe = map(int, f[1:4])
            lo, hi, matches, block, mq = map(int, f[7:12])
            if qs != 0 or qe != qlen or matches != qlen or block != qlen or hi-lo != qlen:
                continue
            if "cg:Z:%d=" % qlen not in f[12:]:
                continue
            key = (f[5], lo, hi, f[4])
            hits[f[0]][key] = max(mq, hits[f[0]].get(key, 0))
    return hits


def locate(j, hits, min_mapq):
    if j["status"] != "awaiting_alignment":
        return j["status"], None
    pair = []
    for side in ("L", "R"):
        found = hits.get(j["join_id"]+"_"+side, {})
        if len(found) != 1:
            return ("unresolved_missing_or_changed_flank" if not found
                    else "unresolved_duplicate_flank"), None
        hit, quality = next(iter(found.items()))
        if quality < min_mapq or quality == 255:
            return "unresolved_mapping_quality", None
        pair.append(hit)
    left, right = pair
    if left[0] != right[0] or left[3] != right[3]:
        return "unresolved_split_or_reoriented_flanks", None
    lo, hi = (left[2], right[1]) if left[3] == "+" else (right[2], left[1])
    if lo >= hi:
        return "unresolved_gap_removed_or_flanks_reordered", None
    return "mapped_pending_gap_check", (left[0], lo, hi, left[3])


def union_bp(intervals):
    end, total = -1, 0
    for lo, hi in sorted(intervals):
        total += max(0, hi-max(lo, end))
        end = max(end, hi)
    return total


def chromosome(alns, lo, hi, minimum, margin):
    per = defaultdict(list)
    for start, end, chrom in alns:
        a, b = max(lo, start), min(hi, end)
        if a < b:
            per[chrom].append((a, b))
    ranked = sorted(((union_bp(iv), ch) for ch, iv in per.items()), reverse=True)
    if not ranked or ranked[0][0] < minimum:
        return "."
    if len(ranked) > 1 and ranked[0][0] < margin*ranked[1][0]:
        return "."
    return ranked[0][1]



def assess(a):
    data = json.loads(Path(a.prefix+".source.json").read_text())
    hits = exact_hits(a.flank_paf)
    refmap = {}
    if Path(a.ref_map).name != "NO_HARMONIZE":
        for r in rows(a.ref_map)[1]:
            name = r["new_name"]
            if name.startswith("chr") and "+" not in name:
                refmap[r["old_name"]] = name.rsplit("_", 1)[0]
    alignments = defaultdict(list)
    if Path(a.ref_paf).name != "NO_PAF":
        opener = gzip.open if a.ref_paf.endswith(".gz") else open
        with opener(a.ref_paf, "rt") as handle:
            for line in handle:
                f = line.rstrip().split("\t")
                chrom = refmap.get(f[5].split("#")[-1]) if len(f) >= 12 else None
                if chrom and int(f[9]) >= a.min_block:
                    alignments[f[0]].append((int(f[2]), int(f[3]), chrom))
    candidates = {r["scaffold"]: r for r in rows(a.candidates)[1] if r["assembly"] == a.assembly}
    by_target = defaultdict(list)
    for j in data["joins"]:
        j["status"], location = locate(j, hits, a.min_mapq)
        j.update(current_object="", current_lo="", current_hi="", orientation="",
                 left_chrom=".", right_chrom=".", junction_assessment="unresolved",
                 auto_cut="no", evidence_profile="not_selected")
        if location:
            name, lo, hi, orient = location
            j.update(current_object=name, current_lo=lo, current_hi=hi, orientation=orient)
            by_target[name].append(j)
    digest = checksum(a.current)
    fields, called = read_table(a.native)
    fields += [k for k in ("transition_lo", "transition_hi", "evidence_only", "chromosome_member", "location_status") if k not in fields]
    if "assembly_sha256" not in fields:
        raise ValueError("Native table must carry current FASTA provenance")
    if any(r["assembly_sha256"] != digest for r in called):
        raise ValueError("Native calls belong to a different FASTA")
    transition_rows = []
    seen = set()
    for name, seq in records(a.current):
        if name in seen:
            raise ValueError("Duplicate current FASTA ID")
        seen.add(name)
        for j in by_target.get(name, []):
            lo, hi = j["current_lo"], j["current_hi"]
            if hi > len(seq) or hi-lo != j["source_hi"]-j["source_lo"] or set(seq[lo:hi]) != {"N"}:
                j["status"] = "unresolved_gap_changed_or_replaced"
                continue
            j["status"] = "recovered_exact_flanks_and_gap"
            alns = alignments.get(name, [])
            lc = chromosome(alns, max(0, lo-a.window), lo, a.min_bp, a.margin)
            rc = chromosome(alns, hi, min(len(seq), hi+a.window), a.min_bp, a.margin)
            j.update(left_chrom=lc, right_chrom=rc,
                     junction_assessment=("unassigned_reference_flank" if "." in (lc, rc) else
                                          "different_reference_chromosomes" if lc != rc else
                                          "same_reference_chromosome"))
    # Enrich only intervals already located on inferred chromosome candidates.
    # Exact recovery audits all old gaps but never creates an independent candidate.
    for ti, row in enumerate(called, 1):
        meta = candidates.get(row["scaffold"], {})
        if meta.get("chromosome_member") != "yes" or "+" not in meta.get("name", ""):
            raise ValueError("Native call outside inferred chromosome candidate scope")
        if row.get("transition_lo", ".") == ".":
            continue
        lo, hi = int(row["transition_lo"]), int(row["transition_hi"])
        matched = [j for j in by_target.get(row["scaffold"], [])
                   if j["status"] == "recovered_exact_flanks_and_gap"
                   and lo <= (j["current_lo"]+j["current_hi"])//2 <= hi]
        safe = [j for j in matched
                if (j["left_chrom"], j["right_chrom"]) ==
                   (row["left_chrom"], row["right_chrom"])]
        native_pos = int(row["cut_bp"]) if row["callable"] == "yes" else None
        # Native AGP midpoints use a 1-based formula; recovered intervals are
        # 0-based. An odd-sized gap can otherwise appear twice one base apart.
        positions = {(native_pos if native_pos is not None and j["current_lo"] < native_pos < j["current_hi"]
                      else (j["current_lo"]+j["current_hi"])//2) for j in matched}
        if native_pos is not None:
            positions.add(native_pos)
        if len(positions) > 1 or (native_pos is not None and matched and not safe):
            row.update(callable="no", candidate_verdict="REVIEW", evidence_only="yes",
                       location_status="unresolved", reason="multiple_current_or_recovered_gaps")
        elif (native_pos is None and len(safe) == 1 and len(positions) == 1
              and row.get("reason") == "no_gap_between_alignment_anchors"
              and hi-lo <= a.transition_max_bridge):
            j = safe[0]
            pos = (j["current_lo"]+j["current_hi"])//2
            row.update(cut_bp=str(pos), agp_join_bp=str(pos), agp_join_distance="0",
                       agp_source="round1_recovered_by_exact_flanks", gap_len=str(j["current_hi"]-j["current_lo"]),
                       callable="yes", candidate_verdict="REVIEW", evidence_only="no",
                       location_status="unique_supported_gap",
                       reason="exact_recovered_gap_between_alignment_anchors;review_required",
                       join_scope="recovered_round1_review_only")
            j["evidence_profile"] = "selected_review"
        transition_rows.append(dict(
            assembly=a.assembly, scaffold=row["scaffold"], transition_id=f'{row["scaffold"]}:T{ti}',
            transition_lo=lo, transition_hi=hi, left_chrom=row["left_chrom"], right_chrom=row["right_chrom"],
            profile_allowed=True, recovered_join_ids=";".join(j["join_id"] for j in matched) or ".",
            recovered_join_count=len(matched), diagnostic_bp=int(row["cut_bp"]),
            profile_status="selected_review", auto_cut="no", assembly_sha256=digest))

    if set(by_target)-seen:
        raise ValueError("Flank alignment target absent from current FASTA")
    write_called(a.prefix+".review_joins.tsv", called)
    columns = ["join_id", "source_object", "source_lo", "source_hi", "status",
               "current_object", "current_lo", "current_hi", "orientation",
               "left_chrom", "right_chrom", "junction_assessment", "auto_cut", "evidence_profile",
               "source_sha256", "current_sha256"]
    with open(a.prefix+".older_join_audit.tsv", "w") as out:
        writer = csv.DictWriter(out, fieldnames=columns, extrasaction="ignore",
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for j in data["joins"]:
            writer.writerow(dict(j, source_sha256=data["source_sha256"], current_sha256=digest))
    # Separate precise location recovery from biological adjudication. This table
    # deliberately does not pretend that read evidence has already been assessed.
    for decision, row in zip(transition_rows, [r for r in called if r.get("transition_lo", ".") != "."]):
        decision.update(
            coordinate_status=row.get("location_status", "unresolved"),
            coordinate_reason=row.get("reason", "."),
            cut_location_eligible=row.get("callable", "no"),
            structural_status=row.get("structural_status", "alignment_transition_requires_adjudication"),
            read_status="not_assessed_in_join_recovery",
            action="REVIEW", action_reason="junction_evidence_adjudication_pending",
            gap_origin=row.get("agp_source", "."),
            decision_scope=row.get("decision_scope", "."))
        for key in ("left_arm_aligned_bp", "right_arm_aligned_bp", "left_arm_fraction",
                    "right_arm_fraction", "assigned_union_bp", "arm_count", "arm_pattern"):
            decision[key] = row.get(key, ".")
    transition_fields = ["assembly", "scaffold", "transition_id", "transition_lo", "transition_hi",
                         "left_chrom", "right_chrom", "profile_allowed", "recovered_join_ids",
                         "recovered_join_count", "diagnostic_bp", "profile_status", "auto_cut",
                         "assembly_sha256", "coordinate_status", "coordinate_reason",
                         "cut_location_eligible", "structural_status", "read_status", "action",
                         "action_reason", "gap_origin", "decision_scope", "left_arm_aligned_bp",
                         "right_arm_aligned_bp", "left_arm_fraction", "right_arm_fraction",
                         "assigned_union_bp", "arm_count", "arm_pattern"]
    with open(a.prefix+".transition_intervals.tsv", "w") as out:
        writer = csv.DictWriter(out, fieldnames=transition_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(transition_rows)
    counts = defaultdict(int)
    for j in data["joins"]:
        counts[j["status"]] += 1
    Path(a.prefix+".older_join_summary.json").write_text(json.dumps(dict(
        assembly=a.assembly, total=len(data["joins"]), statuses=dict(counts),
        transition_intervals=len(transition_rows),
        diagnostic_profiles=sum(r["profile_status"] == "selected_review" for r in transition_rows),
        auto_cut=False, coordinate_system="0-based half-open",
        source_sha256=data["source_sha256"], current_sha256=digest), indent=2)+"\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["prepare", "assess"])
    p.add_argument("--prefix", required=True)
    p.add_argument("--round1-agp")
    p.add_argument("--round1-fasta")
    p.add_argument("--flank", type=int, default=10000)
    p.add_argument("--minimum", type=int, default=1000)
    p.add_argument("--assembly")
    p.add_argument("--current")
    p.add_argument("--flank-paf")
    p.add_argument("--ref-paf")
    p.add_argument("--ref-map")
    p.add_argument("--candidates")
    p.add_argument("--native")
    p.add_argument("--min-mapq", type=int, default=30)
    p.add_argument("--window", type=int, default=500000)
    p.add_argument("--min-bp", type=int, default=100000)
    p.add_argument("--margin", type=float, default=2.0)
    p.add_argument("--min-block", type=int, default=2000)
    p.add_argument("--min-span", type=int, default=0,
                   help="deprecated and ignored; scope comes from candidate chromosome membership")
    p.add_argument("--transition-max-bridge", type=int, default=5000000)
    a = p.parse_args()
    if a.transition_max_bridge < 0 or a.flank < a.minimum or a.minimum < 1 or a.min_bp < 1 or a.window < a.min_bp or a.margin <= 1:
        p.error("Invalid flank/evidence thresholds")
    if a.mode == "prepare":
        prepare(a.round1_agp, a.round1_fasta, a.flank, a.minimum, a.prefix)
    else:
        assess(a)


if __name__ == "__main__":
    main()
