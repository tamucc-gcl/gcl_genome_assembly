# Full-cohort R comparison acceptance — 2026-10-28

The optional finalized-assembly comparison subgate is accepted. No further
comparison rerun or production code fix is required by this packet. This closes
one category 7 subgate, not all operational cleanup or category 6 reporting.

## Execution and inputs

The supplied archive exits 0: 13 completed tasks, zero failed/ignored/retried.
Five mapping jobs ran concurrently. Total elapsed workflow time was about seven
minutes. Mapping task realtime ranged from 1m44s to 6m01s; peak RSS 15.1–33.7 GB.
R plotting took 2.8–9 seconds per pair. First-launch resume and unmatched selector
warnings are harmless for the bounded comparison entry.

Ten final assemblies form exactly five within-sample, same-species pairs. All
ten packaged FAI fingerprints, lengths, sequence counts and source FASTA paths
reconcile with the manifest. These dictionaries are byte-identical to those in
the earlier sparse comparison packet; no assembly change explains the recovered
synteny. Whole-FASTA checksums were not independently recomputed from this packet.

All recorded commands use minimap2 2.30-r1287 asm5 without base-level `-c`, with
a single reference index. All plots use native R 4.5.3, ggplot2 4.0.1, pafr 0.0.2
and gggenomes 1.2.0. All ten figure statuses are `plotted`, and the figures were
visually inspected, including the fragmented CPla and composite CTlk cases.

## Mapping coverage

Fractions below use the complete final assembly denominators. Retained means
MAPQ >=5 and mapping-chain block length >=150 kb.

| Pair | Unfiltered hap1 span coverage | Unfiltered hap2 span coverage | Retained hap1 span coverage | Retained hap2 span coverage |
|---|---:|---:|---:|---:|
| CBau | 97.59% | 99.64% | 96.43% | 96.27% |
| CLim | 96.70% | 98.99% | 94.07% | 90.47% |
| CMat | 96.80% | 99.72% | 95.88% | 96.44% |
| CPla | 87.81% | 98.12% | 75.00% | 81.63% |
| CTlk | 98.02% | 98.42% | 95.81% | 86.68% |

The lower retained coverage is explicitly explained in part by filtering, rather
than disappearing measurements. With MAPQ >=5 and a 10 kb minimum, CPla retains
87.37%/96.78% and CTlk 97.87%/95.86%. Existing 150 kb overview defaults remain
unchanged in this acceptance; the raw mappings and threshold tables preserve
the evidence needed for a later display choice.

All reported fractions were independently checked against covered-base totals
and denominators. Bounds, monotonic filtering and QC-summary totals reconcile.
Raw PAFs are deliberately omitted from the archive, so interval unions cannot
be independently recomputed from this packet. The union algorithm has separate
regression tests. Chain-span coverage is not base-by-base identity or validation
of every adjacency within a long mapping chain.

## Figures and interpretation limits

CBau/CLim/CMat show broad chromosome correspondence. CTlk displays its split and
composite sequence structure and a reversed composite mapping; these are evidence
for biological review, not a reason to reopen accepted manual cuts in this
pipeline acceptance pass. CPla's fragmented assemblies give a scattered dotplot
and dense ribbon view; whole-scaffold order is not equivalent to chromosome
order, so its apparent crossings are not automatically rearrangements.

Plots correctly apply the additional 1 Mb sequence display threshold. The numbers
drawn are 218/305/180/484/289 for CBau/CLim/CMat/CPla/CTlk, versus retained QC
counts 225/327/185/1282/325. These differences are display exclusions, not lost
records. Plotted sequence lengths independently match the final indexes above
the threshold.

Dense/small-scaffold labels, whitespace and subtitle wording remain reasonable
category 6 presentation refinements. The current output is functional and the
sparse-alignment regression is resolved. The approximate match fraction column
must not be interpreted as nucleotide identity; it is explicitly marked as an
approximate mapping-chain measurement.

## Remaining work

Continue category 7's broader dependency/provenance review and final combined
integration, followed by category 6 unified Markdown reporting and scheduled R
figure migration. Categories 1–5 remain accepted. Do not request another identical
pairwise run without a concrete change or newly identified failure.
