# Production graph acceptance — 2026-10-25

Packet: `core-graph-production-02.review.tar.gz`.
Decision: accept the core production graph policy and close category 5. Keep
pre-existing low-retention flags visible for subsequent graph QC and interpretation;
acceptance does not certify complete representation of every assembly base.

## Execution and production policy

Exit 0; 17 completed tasks, 840 cached tasks, no failures, retries or ignored errors.
Cactus construction took 16h 16m elapsed. Source-coordinate audit took 8m 24s;
graph statistics 20m 10s, sharing 32m 49s, similarity 16m 15s and regional
attribution 12m 50s. These downstream tasks overlap where their dependencies permit.
All ten assemblies and the retained CMat hap1 reference are unchanged inputs.

The native Crest extraction/header check passed. One guarded deletion-filter pass
and one corrected regional split ran with `-u 25000 -s`, the existing coverage,
uniqueness and MAPQ settings, and explicit reference chromosomes. The adapter
normalized 12,273,232 mapping records covering 987 distinct fragment identifiers.
No source-coordinate or alias-collision failure occurred. The production path ran
without the temporary assignment/filter/gate/chunk replay branches.

## Recovered chromosome arms

| CTlk source | Input bp | Previous CLIP path bp | New CLIP path bp | New/input |
|---|---:|---:|---:|---:|
| hap1 chr7_1+chr12_1 | 75,393,929 | 40,638,091 | 72,996,680 | 96.82% |
| hap2 chr7_1+chr12_1 | 83,428,940 | 40,814,826 | 79,994,370 | 95.88% |

The production chr7 partitions match the selected diagnostic exactly:
hap1 [0,32,828,442), hap2 [43,248,070,83,428,940). The adjacent chr12 partitions
cover hap1 [32,828,442,75,343,827) and hap2 [40,321,43,248,070). All assigned
intervals survive into FULL. The only unassigned ends on these composites are
50,102 bp in hap1 and 40,321 bp in hap2, matching the diagnostic expectation.
These are graph partitions, not additional assembly cuts.

Separate CTlk chr4 and chr14 paths remain represented. All 134 source sequences
classified as chromosome-scale have a retained CLIP path. No chromosome-scale
source gains a new below-80% retention flag. The largest CLIP path decrease on an
ordinary chromosome relative to the v8 graph is 103,869 bp (CMat hap2 chr6_1,
0.137% of its 76.06 Mb input). Both fragmented CPla assemblies gain about 57 Mb
of aggregate retained path length. These are representation changes, not new
biological sequence or a claim about assembly completeness.

## Source accounting and graph checks

All 6,182 sources reconcile to 10,304,718,010 input bases:

| Preprocessing category | bp |
|---|---:|
| Assigned to chromosome partitions | 9,847,459,478 |
| Ambiguous, excluded from chromosome subproblems | 203,087,437 |
| Outside reported partitions, including softclipped ends | 254,171,095 |

The latter two categories total 457,258,532 bp (4.44% of input). They are explicit
preprocessing exclusions, remain in the source assemblies and are not silently
included in an "other" graph bin. The earlier diagnostic summarized a selected
1,602-source subset; the production totals cover all 6,182 sources and must not be
compared as though their denominators were identical.

Assigned source spans missing from FULL: **0 bp**. FULL spans outside assigned
partitions: **0 bp**. CLIP spans outside FULL: **0 bp**. Original-source intervals
are in bounds and non-overlapping. Distinct graph segment IDs are unique.

| Graph | Distinct nodes | Distinct graph bp |
|---|---:|---:|
| FULL | 146,822,171 | 2,372,383,979 |
| CLIP | 146,033,175 | 1,959,888,337 |

Distinct-node totals reconcile with the existing graph/clipping reports. Sharing
and regional histogram/attribution checks pass for ten haplotypes/five individuals;
variant export/identity checks pass (32,385,711 records, six VCF columns representing
five individuals). The similarity report has zero negative inertia and zero
negative NJ edges. These checks do not certify biological variant interpretations.

## Remaining flags and scope

Five low-retention chromosome paths remain in `source_report.md`, down from seven
when the same screen is applied to the previous v8 graph. All five predate this change:

| Assembly | Source | Previous CLIP bp | Current CLIP bp |
|---|---|---:|---:|
| CBau hap1 | chr8_1 | 61,337,691 | 61,320,976 |
| CLim hap2 | chr8_1 | 59,751,574 | 59,717,872 |
| CMat hap2 | chr8_1 | 59,696,260 | 59,700,471 |
| CTlk hap2 | chr8_1 | 60,811,828 | 60,784,064 |
| CTlk hap2 | chr9_1 | 7,335,077 | 7,335,077 |

The chr8 paths retain about 75–80% of input; the CTlk chr9 piece retains 48.90%.
Their causes are not established by this packet. The source audit correctly remains
REVIEW for them. They are not evidence that the new arm-recovery policy failed,
and do not justify more assembly cuts. Preserve these flags for later graph
self-mapping/representation QC and biological interpretation.

The regional figures correctly retain chr7/chr12 in their own composite row;
grey CTlk chr12 cells marked C describe the naming group, not missing sequence.
Report consolidation/readability, operational resource/output cleanup and removal
of default-off temporary replay/trace diagnostics remain categories 6/7. Keep the
production adapter, pinned corrected tool and permanent source/node audits.
No further graph rebuild is required for acceptance of this change.
