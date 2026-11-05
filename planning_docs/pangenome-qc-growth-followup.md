# Next pangenome work: graph QC, MultiQC and growth fitting

Agreed scope, 2026-11-05. This schedules additions for the next pangenome work;
it does not enable them or reopen the accepted core gate. The canonical checklist
remains [remaining work](remaining-pipeline-work-2026-10-15.md).

## Graph-intrinsic QC and pangenome MultiQC

Implement the measurements missing from the current core, alongside existing
size/count, representation/clipping and variant audits:

- Explicit whole-graph/component acyclicity and self-loop observations. Cycles
  can represent real inversions and must not automatically fail graph QC.
- Node degree distributions and mean/maximum summaries.
- A properly scoped graph node-depth/path-coverage summary, with definitions;
  do not reproduce the old first-20-lines diagnostic as a genome-wide result.
- Mean link length and summed path-node distances, with layout/order provenance.
  Compute on disposable optimized copies without modifying the accepted graph.
- Nucleotide composition and the additional whole-genome/per-chromosome odgi
  statistics exported by the old MultiQC producer.
- Cactus assembly-to-graph GAF match/block-length edit rate and identity,
  aligned bases and alignment counts. Validate field availability and denominator,
  flag secondary/overlapping alignments and distinguish this from independent
  read-to-graph self-mapping. Do not treat summed alignment bp as unique coverage.
- A **pangenome MultiQC HTML report plus its data directory**, incorporating
  graph statistics and current `bcftools stats` inputs. Link it from the Markdown
  pangenome report and main evidence index; it supplements those reports.

Keep species, graph role (FULL/CLIP), chromosome, input identity and method/tool
version explicit. Reuse existing metrics instead of computing duplicate versions.
Publish structured measurements and raw outputs. Missing/unsupported measurements
get explicit statuses; they must not quietly disappear into successful tasks.
Pin supported tools/MultiQC dependencies and test parser compatibility. Validate
small known graphs and a bounded native run before accepting the extension.

These requirements extend categories 10/13. Intrinsic graph QC does not require
read-to-graph mapping, and the two must remain independently schedulable.

## Growth-model fit and open/closed summary

Reimplement the fitted metric for both biological groupings already supported by
Panacus: haplotypes and individuals. Individual membership is the union of its
contributing haplotype paths, counted once per individual. Record haplotype
availability so incomplete individuals are visible. Never pool species or count
synthetic reference paths as biological samples.

Use the descriptive cumulative power-law model `P(m) = k * m^gamma`, preserving
explicit model identity. Fit the union accumulation curve from retained Panacus
results; do not fit every quorum/core curve under the same openness interpretation.
The initial estimator follows the previous log-log fit, with its fit domain and
residual diagnostics recorded. Compute separately from the R/ggplot renderer.

Planned minimum: **five observed units per grouping**, configurable with a hard
floor of three for this two-parameter fit. Five is an execution/reporting guard,
not evidence that population inference is well powered. Require at least three
usable, distinct positive cohort sizes, finite positive accumulation values and
a full-rank fit. Flat curves are legitimate, not a reason to discard a fit.
Below the minimum or on failure, retain the observed curve and emit a clear
not-fitted reason; do not invent a classification or fail unrelated analyses.

Report `k`, `gamma`, unit, cohort size, number/domain of fitted points, fit quality,
graph role and status. Restore the old operational classification rule initially:
`gamma > 0.1` = open; otherwise closed. Make the cutoff explicit/configurable,
record its value and call this a **model-based classification at the configured
cutoff**, not a universal mathematical boundary or proof of future saturation.
Do not conflate this cumulative exponent with exponents fitted to incremental
novel-sequence yield. No extrapolation is required for the initial implementation.

Add one graph-summary row per grouping with exponent, open/closed label, N and
fit status; link coefficients, observed/fitted curves and diagnostics. Report both
small sample size and the dependence of haplotypes within individuals. Avoid
ordinary regression confidence intervals that pretend the nested mean-curve
points are independent experimental replicates.

## Confidence-band exploration: agreed next-pangenome scope

Exploring confidence/variability bands is explicitly included in the next
pangenome visit, alongside graph QC, MultiQC and growth fitting. Evaluate both
haplotype and individual groupings. Distinguish (1) subset/order variability
conditional on the observed cohort, (2) uncertainty in fitted parameters and
curves, and (3) any justified inference beyond the observed cohort. Compare
appropriate analytical and resampling approaches on known examples and real
outputs; preserve sequence co-occurrence and individual grouping where relevant.
Assess small-cohort behavior, computational cost and interval coverage before
choosing which bands to publish and how to label them. Do not automatically
restore the old band or treat exploration as an already-validated method.

### The old uncertainty band and potential replacement

The old renderer used a coverage histogram `h_i` (bp found in exactly i of N
haplotypes). For a random subset of m haplotypes:

- `p_i = 1 - choose(N-i, m) / choose(N, m)`;
- mean accumulated sequence = `sum(h_i * p_i)`;
- purported SD = `sqrt(sum(h_i * p_i * (1-p_i)))`;
- band = mean plus/minus 1.96 times that SD.

The expectation is compatible with without-replacement subset sampling, but the
variance drops covariance between bases. Bases on the same presence/absence
block enter subsets together; they are not independent trials. The histogram
also loses which haplotypes carry each block, so it cannot recover all those
covariances. This is a flaw in the band's interpretation, not a reason to
abandon the mean curve or growth fitting.

For example, if one whole 1-Mb block is present in half the possible subsets,
its true contribution is either 0 or 1 Mb: SD = 500,000 bp. Treating its million
bases as independent instead gives SD = 500 bp. This illustrates how the old
band can be far too narrow; it is not a measured discrepancy for this dataset.

Preferred next implementation: evaluate actual sequence unions over subsets or
random sample orders using carrier memberships and graph sequence lengths.
Enumerate all subsets where bounded; otherwise use reproducible sampling with
recorded seed/replicate count. Use haplotypes for the haplotype curve and whole
individual groups for the individual curve. Validate that empirical means match
Panacus expectations and that linked-block examples give the correct spread.
Do not require an impractically dense base-by-sample matrix.

Label resulting quantiles **subset/order variability within this observed cohort**,
not confidence in the unsampled population or in the fitted model. Zero spread at
all N observed units is expected for this conditional calculation. Population or
fit-parameter uncertainty is a separate resampling/method question. The fitted
metric need not wait for that separate confidence-interval methodology; omit an
unvalidated band rather than retain the old formula.

References:

- [Panacus implementation and supported groupings](https://github.com/codialab/panacus)
- [vegan accumulation documentation: conditional variance, permutations and correlation limitations](https://search.r-project.org/CRAN/refmans/vegan/html/specaccum.html)
