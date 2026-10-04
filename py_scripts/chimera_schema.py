"""Stable called-join contract, including empty/reference and recovered tables."""
import argparse
import csv

JOIN_FIELDS = [
    "assembly", "scaffold", "name", "cut_bp", "left_chrom", "right_chrom",
    "left_component", "right_component", "n_components", "n_transitions",
    "agp_join_bp", "agp_join_distance", "agp_source", "gap_len", "callable",
    "reason", "span_bp", "vote", "candidate_verdict", "transition_lo",
    "transition_hi", "evidence_only", "chromosome_member", "location_status",
    "left_anchor_start", "left_anchor_end", "right_anchor_start", "right_anchor_end",
    "minimum_anchor_bp", "compatible_gap_count", "left_arm_aligned_bp",
    "right_arm_aligned_bp", "left_arm_fraction", "right_arm_fraction",
    "assigned_union_bp", "arm_count", "arm_pattern", "structural_status",
    "decision_scope",
]
CALLED_FIELDS = JOIN_FIELDS + ["assembly_sha256", "coordinate_stage", "join_scope"]


def read_table(path):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader((line for line in handle
                                if line.strip() and not line.lstrip().startswith("#")),
                               delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"Missing or duplicate TSV header: {path}")
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError(f"TSV row width does not match header: {path}")
    return fields, rows


def write_called(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CALLED_FIELDS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            unknown = set(row) - set(CALLED_FIELDS)
            if unknown:
                raise ValueError(f"Unknown called-join columns: {sorted(unknown)}")
            writer.writerow({key: row.get(key, ".") for key in CALLED_FIELDS})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--empty", required=True, help="Write an empty canonical table")
    write_called(parser.parse_args().empty, [])
