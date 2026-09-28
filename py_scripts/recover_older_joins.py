#!/usr/bin/env python3
"""Conservative older-join recovery using exact, unique minimap2 flank alignments.
Coordinates are 0-based half-open internally. Every input join gets an audit row.
Recovered joins are REVIEW only, never automatic cutting authorization.
"""
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



def transition_intervals(alns, length, window, minimum, margin, max_bridge):
    """Coarse uncertainty intervals between successive informative windows.
    Empty/ambiguous windows are bridged only for reporting; broad intervals
    remain visible but do not trigger profiles. No breakpoint is inferred.
    """
    previous = None
    result = []
    for lo in range(0, length, window):
        hi = min(length, lo + window)
        chrom = chromosome(alns, lo, hi, minimum, margin)
        if chrom == ".":
            continue
        if previous and previous[2] != chrom:
            result.append(dict(transition_lo=previous[0], transition_hi=hi,
                               left_chrom=previous[2], right_chrom=chrom,
                               profile_allowed=(lo - previous[1] <= max_bridge)))
        previous = (lo, hi, chrom)
    return result


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
    fields, called = rows(a.native)
    fields += ["transition_lo", "transition_hi", "evidence_only"]
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
            if "." in (lc, rc) or lc == rc:
                continue
            if len(seq) < a.min_span:
                j["evidence_profile"] = "below_scaffold_size_threshold"
                continue
            j["evidence_profile"] = "selected_review"
            meta = candidates.get(name, {})
            # Separate each junction's evidence location. Scaffold-wide votes are not cut authorization.
            row = {k: "." for k in fields}
            row.update(assembly=a.assembly, scaffold=name, name=meta.get("name", name),
                       cut_bp=str((lo+hi)//2), left_chrom=lc, right_chrom=rc,
                       left_component=j["join_id"]+"_current_left",
                       right_component=j["join_id"]+"_current_right",
                       n_components="2", n_transitions="1", agp_join_bp=str((lo+hi)//2),
                       agp_join_distance="0", agp_source="round1_recovered_by_exact_flanks",
                       gap_len=str(hi-lo), callable="yes",
                       reason="exact_coordinate_recovery;biological_decision_requires_review",
                       span_bp=str(len(seq)), vote="not_recomputed_per_junction",
                       candidate_verdict="REVIEW", assembly_sha256=digest,
                       coordinate_stage="pre_finishing", join_scope="recovered_round1_review_only")
            if not any(r["scaffold"] == name and r["cut_bp"] == row["cut_bp"] for r in called):
                called.append(row)

        # Diagnose chromosome changes independently of whether a safe gap exists.
        for ti, transition in enumerate(transition_intervals(
                alignments.get(name, []), len(seq), a.window, a.min_bp,
                a.margin, a.transition_max_bridge), 1):
            lo, hi = transition["transition_lo"], transition["transition_hi"]
            matched = [j for j in by_target.get(name, [])
                       if j["status"] == "recovered_exact_flanks_and_gap"
                       and lo <= (j["current_lo"] + j["current_hi"]) // 2 < hi]
            selected = transition["profile_allowed"] and len(seq) >= a.min_span
            transition_rows.append(dict(
                assembly=a.assembly, scaffold=name, transition_id=f"{name}:T{ti}",
                **transition, recovered_join_ids=";".join(j["join_id"] for j in matched) or ".",
                recovered_join_count=len(matched), diagnostic_bp=(lo+hi)//2,
                profile_status=("selected_review" if selected else
                                "below_scaffold_size_threshold" if len(seq) < a.min_span else
                                "unassigned_span_exceeds_bridge_limit"),
                auto_cut="no", assembly_sha256=digest))
            if not selected:
                continue
            row = {k: "." for k in fields}
            row.update(assembly=a.assembly, scaffold=name,
                       name=candidates.get(name, {}).get("name", name),
                       cut_bp=str((lo+hi)//2), left_chrom=transition["left_chrom"],
                       right_chrom=transition["right_chrom"], callable="no",
                       reason="transition_interval_midpoint;not_a_breakpoint_or_safe_cut",
                       span_bp=str(len(seq)), vote="not_recomputed_per_junction",
                       candidate_verdict="REVIEW", assembly_sha256=digest,
                       coordinate_stage="pre_finishing", join_scope="transition_diagnostic_only",
                       transition_lo=str(lo), transition_hi=str(hi), evidence_only="yes")
            # If a native/gap profile already occupies this exact coordinate, retain it.
            if not any(r["scaffold"] == name and r["cut_bp"] == row["cut_bp"] for r in called):
                called.append(row)

    if set(by_target)-seen:
        raise ValueError("Flank alignment target absent from current FASTA")
    with open(a.prefix+".review_joins.tsv", "w") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(called)
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
    transition_fields = ["assembly", "scaffold", "transition_id", "transition_lo", "transition_hi",
                         "left_chrom", "right_chrom", "profile_allowed", "recovered_join_ids",
                         "recovered_join_count", "diagnostic_bp", "profile_status", "auto_cut",
                         "assembly_sha256"]
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
    p.add_argument("--min-span", type=int, default=20000000)
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
