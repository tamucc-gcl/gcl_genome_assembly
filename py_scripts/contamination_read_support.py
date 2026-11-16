#!/usr/bin/env python3
"""Measure indexed read support around contamination-review intervals."""

import argparse
from collections import defaultdict
from itertools import accumulate
from pathlib import Path

from contamination_evidence import add_read_support_guide, table, write_table


MIN_MAPQ = 20
MIN_FLANK = 50
MAX_GAP = 20
EXCLUDED_FLAGS = 0xF04  # unmapped, secondary, QC failure, duplicate, supplementary


def spans_boundary(blocks, boundary, flank=MIN_FLANK):
    """Require aligned bases on both sides, allowing short indels within the flanks."""
    if not blocks:
        return False
    left, right = boundary - flank, boundary + flank
    previous_end = left
    left_bases = right_bases = 0
    for start, end in blocks:
        if end <= left or start >= right:
            continue
        if start - previous_end > MAX_GAP:
            left_bases = right_bases = 0
        left_bases += max(0, min(end, boundary) - max(start, left))
        right_bases += max(0, min(end, right) - max(start, boundary))
        previous_end = end
        if left_bases >= flank - MAX_GAP and right_bases >= flank - MAX_GAP:
            return True
    return False


def measure(bam, locus):
    sequence = locus["sequence"]
    start, end, length = (int(locus[key]) for key in ("interval_start0", "interval_end0", "contig_bp"))
    if bam.get_reference_length(sequence) != length:
        raise ValueError(f"Read alignment length differs from the review FASTA for {sequence}")
    flank = max(1000, round(.12 * (end - start)))
    lo, hi = max(0, start - flank), min(length, end + flank)
    depth_events = [0] * (hi - lo + 1)
    boundaries = [dict(edge=label, boundary0=position,
                       assessable=MIN_FLANK <= position <= length - MIN_FLANK,
                       reads=set())
                  for label, position in (("Left", start), ("Right", end))]
    for read in bam.fetch(sequence, lo, hi):
        if read.flag & EXCLUDED_FLAGS or read.mapping_quality < MIN_MAPQ:
            continue
        blocks = read.get_blocks()
        for block_start, block_end in blocks:
            a, b = max(lo, block_start), min(hi, block_end)
            if a < b:
                depth_events[a - lo] += 1
                depth_events[b - lo] -= 1
        for edge in boundaries:
            if edge["assessable"] and spans_boundary(blocks, edge["boundary0"]):
                edge["reads"].add(read.query_name)

    depth = list(accumulate(depth_events[:-1]))
    bin_bp = max(1, (len(depth) + 499) // 500)
    identity = dict(action_id=locus["action_id"], sequence=sequence)
    bins = (dict(identity, start0=lo + x, end0=lo + min(x + bin_bp, len(depth)),
                 mean_depth=sum(depth[x:x + bin_bp]) / len(depth[x:x + bin_bp]))
            for x in range(0, len(depth), bin_bp))
    counts = (dict(identity, edge=edge["edge"], boundary0=edge["boundary0"],
                   spanning_reads=len(edge["reads"]) if edge["assessable"] else "NA")
              for edge in boundaries)
    return bins, counts


def main():
    import pysam

    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("bam", "index", "review"):
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    loci = {row["action_id"]: row for row in table(args.review / "locus_plot_hits.tsv")}
    depth, boundaries = [], []
    with pysam.AlignmentFile(args.bam, "rb", index_filename=str(args.index), require_index=True) as bam:
        for locus in loci.values():
            bins, counts = measure(bam, locus)
            depth.extend(bins)
            boundaries.extend(counts)
    write_table(args.review / "locus_read_depth.tsv", depth,
                ["action_id", "sequence", "start0", "end0", "mean_depth"])
    write_table(args.review / "locus_read_boundaries.tsv", boundaries,
                ["action_id", "sequence", "edge", "boundary0", "spanning_reads"])
    if boundaries:
        add_read_support_guide(args.review)


if __name__ == "__main__":
    main()
