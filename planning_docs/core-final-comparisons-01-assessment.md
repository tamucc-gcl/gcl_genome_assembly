# Native final comparisons 01 — review, 2026-10-28

Execution and finalized-input routing pass. Comparison completeness and presentation
need one focused follow-up before the optional visualization gate is finalized.
No assembly or core graph rebuild is indicated.

## Verified

The supplied review archive exits 0. All 13 tasks completed with no failures,
ignored tasks or retries. Five minimap2 alignments ran concurrently, followed by
five plot tasks, the plan, collector and acceptance summary. Alignment runtimes
were 4m50s (CPla), 16m29s (CTlk), 18m35s (CLim), 23m09s (CBau), and 23m32s
(CMat). Peak alignment RSS was 18.0–43.3 GB. These are observed values, not a
general resource ceiling for arbitrary genomes.

Ten finalized assemblies form exactly five within-sample haplotype pairs, all
taxid 373251. Packaged final indexes independently match all ten recorded SHA256
fingerprints, sequence counts and base totals. Input FASTA paths match the source
manifest. The archive does not contain whole FASTAs, so this is index/path
reconciliation rather than independent sequence checksum verification.

Native versions are minimap2 2.30-r1287, Python 3.11.18 and Matplotlib 3.9.4.
All ten figures have `plotted` statuses and were visually inspected. Plot counts
below QC totals are explained by the 1 Mb final-sequence plotting threshold.
Missing dv values remain NA rather than invented identity values. Warnings are
limited to first-launch resume and selectors absent from this bounded workflow.

## Remaining concerns

| Sample | Retained records | Sum of alignment block lengths | Block-length sum / reference length |
|---|---:|---:|---:|
| CBau | 336 | 134,693,081 | 12.53% |
| CLim | 260 | 101,029,917 | 10.20% |
| CMat | 321 | 126,430,477 | 11.92% |
| CPla | 295 | 111,305,172 | 10.71% |
| CTlk | 257 | 88,182,294 | 9.25% |

These percentages are block-length sums divided by reference lengths, **not union
coverage**. They expose sparse retained evidence. All PAF records shorter than
150 kb are discarded upstream. The packet has neither prefilter measurements nor
PAFs, so it cannot establish whether the sparse evidence comes mainly from this
filter, the alignment preset, or sequence characteristics. The aligned-block
match fractions of 87.2–89.3% likewise cannot be interpreted as whole-genome
identity or explained conclusively from this packet.

Dotplot titles misleadingly say alignments >= 0 bp because that is the additional
plot-level filter, although the upstream PAF already excludes blocks <150 kb.
Dotplots lack sequence-name axis labels. CPla's hundreds of sequence boundaries
dominate its dotplot and its ribbon overview is difficult to interpret. Ribbon
colors repeat among many scaffolds and should not imply shared chromosome
identity merely because final scaffold names happen to match.

## Bounded follow-up

Measure raw versus retained alignment counts, block lengths and non-overlapping
reference/query coverage at several length thresholds. Preserve enough alignment
evidence to regenerate plots without rerunning minimap2 just to change plot
thresholds. Separate measurement retention from visual filtering and show the
effective filters accurately. If shorter blocks recover useful synteny, adjust
the display defaults from those measurements; otherwise test a more permissive
assembly preset on a representative pair. Do not lower thresholds and declare
coverage solved without the measurement.

Improve chromosome labels and suppress excessive boundary lines for fragmented
assemblies. Full chromosome overview and fragmented-scaffold displays should have
different presentation density while preserving underlying evidence. Keep these
rules based on input structure, not sample names.

Category 7 stays open for this comparison follow-up, remaining dependency audit
and combined integration. Categories 2/3/5 remain accepted; unified reporting
remains category 6.
