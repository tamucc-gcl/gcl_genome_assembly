#!/usr/bin/env python3
"""Restrict extension to explicitly inferred chromosome scaffolds and audit pass-through."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def fingerprints(path):
    records = {}
    name = None
    digest = None
    length = 0
    with open(path, encoding="utf-8") as source:
        for line in source:
            if line.startswith(">"):
                if name is not None:
                    records[name] = dict(length=length, sha256=digest.hexdigest())
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError("Duplicate FASTA ID: " + name)
                digest = hashlib.sha256()
                length = 0
            elif line.strip():
                if name is None:
                    raise ValueError("Sequence before FASTA header")
                bases = line.strip().encode("ascii")
                digest.update(bases)
                length += len(bases)
    if name is not None:
        records[name] = dict(length=length, sha256=digest.hexdigest())
    if not records:
        raise ValueError("Empty assembly FASTA")
    return records


def write_audit(prefix, rows):
    fields = ["scaffold", "original_length", "chromosome_member", "eligible", "reason",
              "final_length", "added_bp"]
    with open(str(prefix) + ".teloclip_scope.tsv", "w", encoding="utf-8", newline="") as dest:
        writer = csv.DictWriter(dest, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def prepare(fasta, name_map, prefix):
    records = fingerprints(fasta)
    mapped = {}
    if Path(name_map).name != "NO_HARMONIZE":
        with open(name_map, encoding="utf-8") as source:
            reader = csv.DictReader(source, delimiter="\t")
            if not {"old_name", "chromosome_member"}.issubset(reader.fieldnames or []):
                raise ValueError("Telomere scope requires the current harmonized map schema")
            for row in reader:
                name = row["old_name"]
                if name in mapped:
                    raise ValueError("Duplicate name-map ID: " + name)
                if row["chromosome_member"] not in {"yes", "no", "unresolved"}:
                    raise ValueError("Invalid chromosome membership for " + name)
                mapped[name] = row
        if set(mapped) != set(records):
            raise ValueError("Telomere name-map IDs do not match the assembly")
    rows = []
    eligible = []
    with open(str(prefix) + ".eligible.bed", "w", encoding="utf-8") as bed:
        for name, record in records.items():
            row = mapped.get(name, {})
            member = row.get("chromosome_member", "unavailable")
            selected = member == "yes"
            if selected:
                eligible.append(name)
                bed.write(f"{name}\t0\t{record['length']}\n")
            reason = ("in_inferred_chromosome_set" if selected else
                      row.get("chromosome_scope_reason") or
                      ("outside_inferred_chromosome_set" if member == "no" else
                       "chromosome_classification_" + member))
            rows.append(dict(scaffold=name, original_length=record["length"],
                             chromosome_member=member, eligible="yes" if selected else "no",
                             reason=reason, final_length=record["length"], added_bp=0))
    Path(str(prefix) + ".scope.json").write_text(json.dumps(dict(records=records, eligible=eligible, audit=rows)), encoding="utf-8")
    write_audit(prefix, rows)
    return len(eligible)


def restore_order(fasta, names):
    """Copy FASTA byte ranges in source order without loading chromosome sequences.

    Teloclip 0.3.2 writes extended records before copying untouched records.
    Ordering is a pipeline responsibility, not a sequence-validity failure.
    """
    path = Path(fasta)
    blocks = {}
    current = None
    with path.open("rb") as source:
        while True:
            start = source.tell()
            line = source.readline()
            if not line:
                if current is not None:
                    blocks[current][1] = start
                break
            if line.startswith(b">"):
                if current is not None:
                    blocks[current][1] = start
                current = line[1:].split()[0].decode("utf-8")
                blocks[current] = [start, None]
    temporary = path.with_name(path.name + ".ordered.tmp")
    try:
        with path.open("rb") as source, temporary.open("wb") as dest:
            for name in names:
                start, end = blocks[name]
                source.seek(start)
                remaining = end - start
                last = b""
                while remaining:
                    chunk = source.read(min(1024 * 1024, remaining))
                    if not chunk:
                        raise ValueError("Truncated FASTA while restoring order")
                    dest.write(chunk)
                    last = chunk[-1:]
                    remaining -= len(chunk)
                if last != b"\n":
                    dest.write(b"\n")
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def verify(fasta, prefix):
    state = json.loads(Path(str(prefix) + ".scope.json").read_text(encoding="utf-8"))
    records = fingerprints(fasta)
    if set(records) != set(state["records"]):
        missing = sorted(set(state["records"]) - set(records))
        unexpected = sorted(set(records) - set(state["records"]))
        raise ValueError(f"Telomere extension changed assembly IDs: missing={missing}, unexpected={unexpected}")
    eligible = set(state["eligible"])
    for row in state["audit"]:
        name = row["scaffold"]
        if name not in eligible and records[name] != state["records"][name]:
            raise ValueError("Telomere extension changed an ineligible sequence: " + name)
        row["final_length"] = records[name]["length"]
        row["added_bp"] = records[name]["length"] - row["original_length"]
    if list(records) != list(state["records"]):
        restore_order(fasta, state["records"])
    write_audit(prefix, state["audit"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--fasta", required=True)
    prep.add_argument("--name-map", required=True)
    prep.add_argument("--prefix", required=True)
    check = sub.add_parser("verify")
    check.add_argument("--fasta", required=True)
    check.add_argument("--prefix", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        print(prepare(args.fasta, args.name_map, args.prefix))
    else:
        verify(args.fasta, args.prefix)


if __name__ == "__main__":
    main()
