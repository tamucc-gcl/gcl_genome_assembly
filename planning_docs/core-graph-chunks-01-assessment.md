# Chunk comparison: accepted diagnostic result and production candidate

Packet `core-graph-chunks-01.review.tar.gz` completed with exit 0 and no Nextflow
warnings/errors. 889 tasks: 875 cached, 14 completed (assembly summary, 12 assignment
jobs and one comparison report). New assignments took 64–68 seconds each; report
took 1.9 seconds. All 12 sets of 1,602 coordinate audits pass (19,224 source checks),
with zero duplicate/overlapping assignments.

## Preferred candidate

Chromosome-isolated deletion filtering, corrected grouping and the existing
`-u 25000` option recover all available post-filter mapped bases of both CTlk
chr7/chr12 composites. This is a supported parameter change plus the previously
isolated grouping correction; no additional upstream decision rule was rewritten.

| Setting, guarded filtering | hap1 retained mapped bp | hap2 retained mapped bp | Interpretation |
|---|---:|---:|---|
| Corrected, 1 Mb | 40,484,757 | 40,602,894 | Neither chr7 arm retained |
| Corrected, 100 kb | 40,484,757 | 79,621,204 | Only hap2 arm recovered |
| Corrected, 25 kb | 72,629,202 | 79,621,204 | Both arms recovered completely relative to filtered input |
| Original grouping, 25 kb | 40,711,058 | 41,773,721 | Only small pieces recovered |

The corrected 25 kb case assigns each recovered chr7 arm as one region:
hap1 [0,32,828,442), hap2 [43,248,070,83,428,940). Their mapped union fractions are
97.92% and 97.11%. These coordinate partitions prepare graph inputs; they are not
new assembly cuts or biological breakpoint decisions.

## Cohort side effects

Across all 1,602 selected sources, the preferred case adds 98,742,583 mapped bases
relative to its 1 Mb baseline. 120 sources gain; none loses mapped coverage.
All 16 ordinary chromosome controls retain mapped coverage and gain no additional
chromosome labels. This is a diagnostic selection rather than an exhaustive test
of every source scaffold.

There are 122 newly assigned source/chromosome combinations:

- Eight include a region at least 1 Mb long: the two CTlk arms and six regions of
  already named CPla composites. Each is assigned as one region.
- 114 have only smaller regions, totaling 18,264,357 mapped bases. Of these, 111
  are previously unassigned sources gaining a chromosome label. Only three add
  a second chromosome label to an already assigned source; all are existing CPla
  composite names. Their added regions have 83.99%, 99.61% and 99.46% mapped union
  coverage. The small-region flag by itself does not establish a false assignment.
- 100 of the 114 small assignments have at least 80% mapped coverage of their
  assigned span; the remaining 14 are at least 50% covered. Their biological
  correctness is not proven by coverage alone.

Native deletion filtering still prevents recovery: at corrected 25 kb, hap1 retains
40,484,757 mapped bases and hap2 43,348,632. Changing only the chunk option cannot
recover mappings already deleted. Disabled filtering also recovers the motivating
arms; the guarded version retains filtering within each chromosome lane.

## Source-end exclusion needs explicit accounting

The chunk option also affects dominant-region smoothing and the existing `-s`
softclip rule. Under guarded/corrected 25 kb, unreported source bases increase from
93,111,000 to 130,916,015: **37,805,015 additional bases** across the selection,
with changes in 112 sources. Ambiguous bases decrease from 240,604,415 to 163,645,845.
Assigned span increases from 1,481,346,544 to 1,520,500,099.

No previously retained mapped coverage is lost, so the extra unreported intervals
do not remove previously surviving mapped bases. Nevertheless, those source bases
exist and cannot be described as fully represented merely because mapped coverage
and overlap audits pass. On the motivating composites, 50,102 hap1 terminal bases
and 40,321 hap2 initial bases are excluded from named partitions; ambiguous span
is zero in both. The current high-level diagnostic table omits unreported span,
so reading its ambiguous column alone would overstate completeness.

## Next production step

Use guarded filtering + corrected grouping + 25 kb as the supported candidate for
the normal graph construction path. Do not run another broad parameter sweep on the
basis of this packet. Integrate the selected preprocessing with exact provenance,
retain the source ownership assertions, then rebuild and check actual FULL/CLIP
source representation and distinct-node accounting. Explicitly account for excluded
source ends and determine how the ordinary unassigned/other route handles them;
do not silently discard them or pad chromosome labels just to hide missing sequence.
Final category acceptance remains open until that real graph run passes. Temporary
replay/trace branches should be removed during operational cleanup after acceptance.
