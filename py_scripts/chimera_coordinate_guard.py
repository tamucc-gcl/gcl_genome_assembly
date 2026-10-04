#!/usr/bin/env python3
"""Validate current AGP/FASTA coordinates and stamp breakpoint-table provenance."""
import argparse
import csv
import hashlib
import re
from pathlib import Path
from chimera_schema import read_table, write_called


def fasta_records(path):
    name, parts = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                if name is None and line.strip():
                    raise ValueError("Sequence before FASTA header")
                parts.append(line.strip())
        if name is not None:
            yield name, "".join(parts).upper()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", required=True)
    p.add_argument("--agp", required=True)
    p.add_argument("--table", required=True)
    p.add_argument("--audit", required=True)
    p.add_argument("--round", required=True)
    p.add_argument("--alignment-status", choices=("assessed", "unavailable"), default="assessed")
    a = p.parse_args()
    lengths, gaps = {}, {}
    with open(a.agp) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            name, lo, hi = f[0], int(f[1]), int(f[2])
            if lo != lengths.get(name, 0) + 1 or hi < lo:
                raise ValueError("Noncontiguous AGP object: " + name)
            lengths[name] = hi
            if f[4] in ("N", "U"):
                if int(f[5]) != hi-lo+1:
                    raise ValueError("AGP gap length mismatch")
                gaps.setdefault(name, []).append((lo-1, hi))
    seen = set()
    for name, seq in fasta_records(a.fasta):
        if name in seen or lengths.get(name) != len(seq):
            raise ValueError("AGP and assessment FASTA differ: " + name)
        seen.add(name)
        for lo, hi in gaps.get(name, []):
            if re.fullmatch("N+", seq[lo:hi]) is None:
                raise ValueError("AGP gap absent in assessment FASTA: " + name)
    if seen != set(lengths):
        raise ValueError("AGP and FASTA sequence sets differ")
    digest = hashlib.sha256()
    with open(a.fasta, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b""):
            digest.update(chunk)
    checksum = digest.hexdigest()
    table = Path(a.table)
    fields, rows = read_table(table)
    if "assembly" not in fields:
        raise ValueError("Breakpoint table has no valid header")
    for row in rows:
        row.update(assembly_sha256=checksum, coordinate_stage="pre_finishing",
                   join_scope=a.round + "_only")
    write_called(table, rows)
    with open(a.audit, "w") as out:
        out.write("metric\tvalue\n")
        out.write("coordinate_stage\tpre_finishing\n")
        out.write(f"assembly_sha256\t{checksum}\n")
        out.write(f"join_scope\t{a.round}_only\n")
        out.write("agp_fasta_check\tpassed\n")
        out.write(f"alignment_transition_assessment\t{a.alignment_status}\n")
        out.write("older_join_provenance\t" +
                  ("unresolved_across_correction;not_tested_for_cuts" if a.round == "round2" else "not_applicable") + "\n")


if __name__ == "__main__":
    main()
