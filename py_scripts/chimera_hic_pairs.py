#!/usr/bin/env python3
"""Lift last-scaffolding-input Hi-C pairs through exactly one matching AGP.

All positions are 1-based. Split components are resolved by source interval.
Ambiguous duplicated placements are discarded, never resolved by arrival order.
This projects existing mappings; it does not realign reads or revise mapping quality.
"""
import argparse
import gzip
from collections import Counter, defaultdict
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pairs", required=True)
    p.add_argument("--round1", required=True, help="AGP for the last scaffolding round")
    p.add_argument("--round2", default="", help="Deprecated: multi-round chains are rejected")
    p.add_argument("--scaffolds", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--pair-type", default="UU")
    a = p.parse_args()
    if a.round2:
        p.error("Use pairs mapped to the last scaffolding input and its single AGP")
    want = set(a.scaffolds.split(","))
    placements, lengths = defaultdict(list), {}
    with open(a.round1) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                raise ValueError("Invalid AGP row: " + line)
            obj, lo, hi = f[0], int(f[1]), int(f[2])
            lengths[obj] = max(lengths.get(obj, 0), hi)
            if f[4] in ("N", "U"):
                continue
            start, end, orient = int(f[6]), int(f[7]), f[8]
            if hi - lo != end - start or orient not in ("+", "-"):
                raise ValueError("AGP component length/orientation is not liftable: " + line)
            placements[f[5]].append((obj, lo, hi, start, end, orient))
    if not want or not want <= lengths.keys():
        raise ValueError("Candidate scaffolds missing from the current AGP")
    dest = Path(a.outdir)
    dest.mkdir(parents=True, exist_ok=True)
    handles = {s: (dest / f"{a.label}.{s}.pairs").open("w") for s in sorted(want)}
    # Keep alternative partners without mixing them into within-scaffold profiles.
    partners = gzip.open(dest / f"{a.label}.alternative_partners.pairs.gz", "wt")
    partners.write("## pairs format v1.0\n#sorted: none\n")
    partners.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n")
    stats = Counter()

    def locate(chrom, pos):
        found = [v for v in placements.get(chrom, []) if v[3] <= pos <= v[4]]
        if len(found) != 1:
            stats["ambiguous_ends" if found else "unmapped_ends"] += 1
            return None
        v = found[0]
        mapped = v[1] + pos - v[3] if v[5] == "+" else v[2] - (pos - v[3])
        if not v[1] <= mapped <= v[2]:
            raise ValueError("Liftover position outside its AGP interval")
        return v, mapped

    opener = gzip.open if a.pairs.endswith(".gz") else open
    try:
        with opener(a.pairs, "rt") as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    raise ValueError("Expected standard pairtools columns")
                stats["pairs_read"] += 1
                if a.pair_type and f[7] != a.pair_type:
                    continue
                x, y = locate(f[1], int(f[2])), locate(f[3], int(f[4]))
                if x is None or y is None:
                    continue
                v1, t1 = x
                v2, t2 = y
                if v1[0] != v2[0]:
                    stats["skipped_cross_scaffold"] += 1
                    if v1[0] in want or v2[0] in want:
                        s1, s2 = f[5], f[6]
                        if v1[5] == "-":
                            s1 = {"+": "-", "-": "+"}.get(s1, s1)
                        if v2[5] == "-":
                            s2 = {"+": "-", "-": "+"}.get(s2, s2)
                        partners.write(f"{f[0]}\t{v1[0]}\t{t1}\t{v2[0]}\t{t2}\t{s1}\t{s2}\t{f[7]}\n")
                        stats["alternative_partner_pairs_written"] += 1
                    continue
                sc = v1[0]
                if sc not in want:
                    continue
                if f[1] == f[3] and v1 == v2:
                    stats["distance_checked"] += 1
                    if abs(t2-t1) != abs(int(f[4])-int(f[2])):
                        raise ValueError("Intra-placement distance changed")
                s1, s2 = f[5], f[6]
                if v1[5] == "-":
                    s1 = {"+": "-", "-": "+"}.get(s1, s1)
                if v2[5] == "-":
                    s2 = {"+": "-", "-": "+"}.get(s2, s2)
                if t1 > t2:
                    t1, t2, s1, s2 = t2, t1, s2, s1
                    stats["swapped"] += 1
                handles[sc].write(f"{f[0]}\t{sc}\t{t1}\t{sc}\t{t2}\t{s1}\t{s2}\t{f[7]}\n")
                stats["pairs_written"] += 1
    finally:
        partners.close()
        for handle in handles.values():
            handle.close()
    for sc in sorted(want):
        (dest / f"{a.label}.{sc}.chromsizes").write_text(f"{sc}\t{lengths[sc]}\n")
    with (dest / f"{a.label}.hic_pairs_audit.tsv").open("w") as out:
        out.write("metric\tvalue\ncoordinate_stage\tpre_finishing\n")
        out.write("transform\tlast_round_AGP_only\n")
        for key in ("pairs_read", "pairs_written", "ambiguous_ends", "unmapped_ends",
                    "skipped_cross_scaffold", "alternative_partner_pairs_written", "distance_checked", "swapped"):
            out.write(f"{key}\t{stats[key]}\n")
        out.write("distance_mismatches\t0\n")
        out.write("status\t" + ("projected" if stats["pairs_written"] else "no_supporting_pairs") + "\n")


if __name__ == "__main__":
    main()
