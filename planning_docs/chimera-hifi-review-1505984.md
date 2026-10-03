# HiFi junction review: run 1505984

Reviewed submitted `chimera-review-20261003-080345-1505984.tar.gz` and
`chimera-hifi-alignments-1505984.tar.gz`. These are the pre-sync baseline, not a
runtime test of the new junction diagnostics. Run succeeded: 11 executed, 289
cached, no failed tasks. SAM records and depth tables were inspected; no pipeline
or Python tests were run locally.

Assessment FASTA SHA256 values recorded by the successful context tasks:

- CTlk hap1: `1fd4707ae4958ce31f7e95b5732e01d70686755c6a62a414b080fac24e4fd96a`
- CTlk hap2: `121ce840edf743fb0a1c9bcc39ef210c241a585bd9e2f6dac3758176bc12ff19`

These match the preceding review's pre-finishing coordinate provenance. The
FASTA files themselves were not transferred or independently rehashed here.

## Midpoint inspection

The reported positions are diagnostic midpoints, NOT established breakpoints.
Count primary mapped records excluding secondary, supplementary, duplicate and
QC-fail flags; require MAPQ >=20 and reference alignment bounds at least 1 kb
on either side of the midpoint. The final column additionally requires
NM / (M + = + X + I + D CIGAR lengths) <=0.01. That is a descriptive screening
criterion for this review, not a validated production threshold or proof of
unique anchors. Insertions/deletions may still occur inside the bracketing record.

| Haplotype | Candidate | Scaffold | Midpoint | Bracketing records | Also <=1% NM rate |
|---|---|---|---:|---:|---:|
| hap1 | 1 | scaffold_1 | 62780042 | 45 | 24 |
| hap1 | 2 | scaffold_1 | 95962670 | 28 | 21 |
| hap1 | 3 | scaffold_1 | 116029144 | 5 | 5 |
| hap1 | 4 | scaffold_1 | 116196428 | 15 | 7 |
| hap1 | 5 | scaffold_5 | 35512640 | 18 | 16 |
| hap2 | 1 | scaffold_5 | 43144255 | 43 | 21 |
| hap2 | 2 | scaffold_7 | 39618934 | 15 | 6 |

No qualifying MAPQ >=20 primary record brackets an entire candidate interval
plus the two 1 kb anchors. Intervals are roughly 50–324 kb, so that absence is
not evidence of an assembly error without evaluating read-length opportunity.
Midpoint support likewise cannot clear the rest of a wide interval.

## Filtered-depth correction

Hap1 scaffold_5 has zero MAPQ >=20 primary depth over
`[35476797,35488017)` (11,220 bp). Seven primary alignments nevertheless bracket
that whole interval plus 1 kb at either end, with MAPQ 1–19. Each aligns
11,198–11,218 of the interval's reference bases as M/= /X operations, with no
overlapping deletion larger than 2 bp. Therefore this is NOT an absence of read
alignments and must not be used as a no-read-support cut. Ambiguous placement
needs investigation; regional SAMs do not include all placements elsewhere in
the assessment assembly.

Hap1 scaffold_1's two 100 bp zero-depth stretches exactly match recovered N gaps
`[116057447,116057547)` and `[116171488,116171588)`. Their recovery audit labels
both flanks chr12. No MAPQ >=20 primary record brackets either gap with 1 kb
flanks. Neither finding establishes the proposed chr12/chr4 transitions.

Some midpoint-bracketing alignments have MAPQ 60 but thousands of NM errors and
supplementary placements. An outer CIGAR span plus MAPQ alone is insufficient
for declaring continuity; preserve CIGAR, NM, clipping and alternative placements.

## Decision and next work

Keep cutting disabled. None of the diagnostic midpoints is justified as a cut by
this inspection. This does not clear the entire composite scaffold.

Use the retained BAMs and candidate FASTAs to locate changes in reliable local
read support throughout each interval, evaluate alignments on both sides of any
proposed boundary, and investigate alternative placements of ambiguous reads.
Do not remap solely to obtain these diagnostics. Candidate FASTAs and all reported
placements for regional read names are the next useful compact export; the full
assessment BAMs should remain on the cluster. Preserve the pre-sync baseline.
