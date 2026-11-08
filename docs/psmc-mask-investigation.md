# Bounded Dipcall mask investigation

This is a removable diagnostic, not a change to production confidence thresholds.
It uses Bash/Slurm rather than another Nextflow run for the alignment audit.
The existing production complete-HiFi route runs independently as a comparator.
There is only one production Nextflow resume, so the two batches can run together.

## Questions measured

- How much reference sequence has no observed alignment span?
- How much is covered only by secondary alignments, scoring-incompatible records,
  alignments shorter than 10 kb, or low-MAPQ alignments?
- Among qualifying records, how much is excluded by overlap? Are the overlapping
  segments from one query sequence or multiple query sequences?
- How much is unique but supported only by a 10-50 kb alignment?
- Does this independently calculated coverage reproduce the actual paftools mask?
- Which haplotype limits the paired mask? Does the hap1 self-alignment already
  lose most sequence, or does the loss arise mainly in hap2 or their intersection?
- Does a complete-HiFi diploid analysis retain more usable sequence, and how does
  its history compare with the assembly route on the same own reference?

The audit follows installed paftools `call` semantics: PAF block length (column 11),
10 kb coverage length, 50 kb variant-support length, MAPQ 5, rejection of secondary
and inversion-secondary records and the missing-s2 scoring check. Categories are
hierarchical and partition selected non-gap scope without double counting.
Unknown/reference gaps and optional input scope masks retain their original
coordinate meaning. Overlap is not automatically labeled a biological repeat.

## Retained versus regenerated evidence

The main pipeline uses scratch. Raw PAFs and per-haplotype masks may have been
removed because they were not published outputs. Preparation checks **only the
exact successful Dipcall work directories in the supplied trace**. Existing PAFs
are reused. Otherwise, the recorded Conda environment runs the standard Dipcall
PAF/mask Make target, with unchanged mapper settings. It does not regenerate
SAM/BAM mappings, VCF calls or assembly-based PSMC fits. These diagnostic outputs
stay in the comparison directory; the original cache/results are read-only.

Reference SHA256, scope agreement, tool/script fingerprints, query/PAF hashes and
commands are retained. The two regenerated masks are compared with the original
published diploid mask. A mismatch is an explicit review failure, not a reason to
quietly substitute the new mask. Missing tasks also produce a failed report and
compact review archive, with available task logs.

## Crest batch

After syncing, run from the project root:

```bash
cd /work/birdlab/GCL/spratelloides_delicatulus_genome
module load miniconda3

bash gcl_genome_assembly/scripts/comparisons/submit_psmc_mask_investigation.sh \
  comparisons/core-psmc-02 \
  comparisons/psmc-mask-investigation-03 \
  Sde-CMat_203

sbatch gcl_genome_assembly/scripts/comparisons/run_core_integration.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-psmc-02/results \
  data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv \
  comparisons/psmc-hifi-comparison-03 \
  --run_report true --report_chromosome_painting true \
  --run_psmc true --psmc_bootstraps 10 \
  --psmc_manifest "$PWD/comparisons/psmc-mask-investigation-03/hifi_manifest.tsv" \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json"
```

The first command prepares ten haplotype tasks for the five-individual baseline,
submits an unthrottled Slurm array (8 CPUs/96 GB per task), then submits a small
collector after all tasks finish, including failures. Slurm controls placement;
there is no additional array concurrency cap. PAF mapping uses the allocated
threads; interval accounting is a sorted endpoint sweep rather than a genome-size
matrix. The report uses the recorded R/ggplot environment, with no runtime installs.
Preparation prints how many PAFs need regeneration. Optional `--analysis-env` and
`--plot-env` prefixes explicitly override recorded environments when needed.

The second command runs only the chosen individual's PSMC HiFi analysis through
the real pipeline. The full own-individual HiFi read set, callability rules,
fitting settings and independent bootstrap jobs are the production implementation.
Scaling remains unset. Core assembly/scaffolding/chimera/graph work should resume
from cache. Mat is an explicit test selection in the command, not hardcoded into
the diagnostic or pipeline. New mapping/calling may take substantially longer
than the mask audit; both batches can proceed in parallel.

Return both:

- `comparisons/psmc-mask-investigation-03.review.tar.gz`
- `comparisons/psmc-hifi-comparison-03.review.tar.gz`

The mask directory contains `README.md`, a stacked exclusion figure, per-contig
tables, classified BED intervals, per-haplotype provenance, and exact paired-mask
counts. Large PAFs remain under `tasks/<individual>.<haplotype>/raw/`, excluded
from the archive. The comparator archive supplies input QC, fitted scope,
heterozygosity and observed/bootstrap histories for comparison with pilot 02.

## Interactive node use

Inside an existing compute-node allocation, after preparation, one task can also
run without Slurm submission or Nextflow:

```bash
bash gcl_genome_assembly/scripts/comparisons/run_psmc_mask_investigation.sbatch \
  "$PWD/comparisons/psmc-mask-investigation-03/manifest.json" task 0
```

Use `input_inventory.tsv` to choose task indices. Once all tasks finish, use the
same script with `report` instead of `task 0` to collect the archive. Do not run
the same task twice concurrently. No cache recovery search is performed.

## Local checks

`tests/test_psmc_mask_investigation.py` tests the filter boundaries, overlapping
records, query identity, gap scope and actual native paftools mask agreement.
`tests/run_psmc_mask_investigation_checks.py` uses a completed native PSMC fixture
to test both PAF reuse and regeneration, Slurm spool-directory script execution,
R rendering, incomplete-task failure and compact archive contents.

## Assembly-preset comparison after the mask audit

Use the completed investigation directory on Crest, where its raw PAFs remain,
not the compact downloaded review archive. After syncing the repository, from
the project root:

```bash
module load miniconda3
bash gcl_genome_assembly/scripts/comparisons/submit_psmc_alignment_presets.sh \
  comparisons/psmc-mask-investigation-03 \
  comparisons/psmc-alignment-presets-04 \
  Sde-CMat_203
```

This submits four unthrottled array tasks: asm10/hap1, asm10/hap2, asm20/hap1,
asm20/hap2, with eight CPUs and 96 GB requested per task. Both asm5 controls
are reused from the successful audit without remapping. The selected individual
is an argument, not a dataset-specific rule. An after-any collector produces
`comparisons/psmc-alignment-presets-04.review.tar.gz`, including failure evidence
if any task fails. The existing HiFi comparison continues independently.

Each trial uses the recorded analysis environment and the generated Dipcall
PAF/mask target. A GNU make command-line override changes only `-xasm5` to the
requested preset in `MM2_OPT`; no installed tool source or production setting
is edited. The 50 kb mask/variant requirement, 10 kb coverage requirement,
MAPQ 5 and overlap rules stay fixed. No SAM/BAM/VCF stage or PSMC fit is run.
The collector verifies matching input hashes, scope, tool versions and filter
rules; predicted masks must still reproduce actual paftools masks exactly.
Gains and losses versus asm5 are measured rather than rejected as discrepancies.

The report includes paired-mask gains/losses, per-haplotype exclusions, block
length counts, block-weighted identity and gap-compressed divergence, commands,
and provenance. Identity summaries include overlapping records and must not be
read as unique genome coverage or heterozygosity. Greater retained coverage is
evidence to examine, not automatic validation of a preset. These diagnostics
produce no new PSMC input: adopting a preset requires rerunning the full Dipcall
alignment/variant route consistently. Merely enlarging a BED mask against an
older VCF could misclassify missed variants as homozygous reference.

`tests/run_psmc_preset_checks.py` runs the four real native preset trials on a
small completed mask fixture, checks actual minimap2 command logs, verifies
unchanged controls and mask reconstruction, renders R/ggplot figures, and checks
incomplete-task status and compact failure archives. It accepts `--source` and
`--out` directories and uses the environments recorded in the source manifest.

## Full variant calling and PSMC comparison

After the successful preset-mask comparison, run:

```bash
module load miniconda3
bash gcl_genome_assembly/scripts/comparisons/submit_psmc_full_presets.sh \
  comparisons/psmc-alignment-presets-04 \
  comparisons/psmc-full-presets-05
```

This prepares the audited individual and scope automatically and submits two
independent Nextflow controllers. Each invokes the **production PSMC workflow**
through the removable `psmc_from_final.nf` comparison entry, using accepted final
eligibility/FASTA inputs. Core assembly/graph work is not invoked. This is a
bounded test of the same modules used by `main.nf`, not a second implementation
of Dipcall, masking or fitting. The normal pipeline exposes the same
`psmc_assembly_preset` parameter.

Both PAF/mask and SAM/VCF targets are regenerated consistently for each preset,
followed by PSMC input filtering, observed fitting, ten parallel-scheduled
bootstrap replicates, export, and R figures/reports. Other fitting settings are
inherited from baseline run provenance, input filters and optional scope BEDs
from the audited source. Scaling remains null. The preparation checks that the
baseline final assembly paths match those tested. Own-reference verification
then runs normally. New VCF/BED outputs remain available outside the archive.

Separate launch/work directories avoid sharing an active Nextflow session with
each other or the HiFi comparison. This intentionally does not resume the main
pipeline session; no expensive upstream assembly jobs are registered. Resources,
environments and per-process Slurm scheduling come from the production config.
Bootstraps are separate jobs under the existing `psmc_bootstrap_max_forks` limit.

Return both archives:

- `comparisons/psmc-full-presets-05/asm10.review.tar.gz`
- `comparisons/psmc-full-presets-05/asm20.review.tar.gz`

Review final callability and exclusion reasons, heterozygosity and curve changes
against asm5 and the independent HiFi comparison before choosing a production
default. More mask coverage alone is insufficient to validate the newly added
sequence. The main default remains unchanged while testing.

## Genotype reliability on shared and preset-specific scope

Use the completed full runs and, optionally, an existing completed HiFi mapping:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_psmc_reliability.sbatch \
  comparisons/psmc-full-presets-05 \
  comparisons/psmc-reliability-06 \
  --hifi-run comparisons/psmc-hifi-comparison-03
```

The HiFi run may still be calling variants; its mapping must have completed.
Only the exact successful mapping work directory recorded in its trace is used.
The own-reference hash is checked against both full assembly runs and the
retained mapping reference. No task is stopped, modified or rerun. Inputs are
read-only; the audit does not require the large files in a downloaded archive.

This removable diagnostic streams the two retained masked diploid consensuses,
preserving exact reference coordinates. It counts shared, asm10-only, asm20-only
and neither-callable regions; reproduces both callable totals; compares actual
unphased genotypes; and reports heterozygosity by region. Shared heterozygote-site
agreement is reported separately so agreeing homozygous-reference calls cannot
hide differences relevant to PSMC. One-Mb window tables and an R/ggplot figure
show where disagreement occurs.

With a HiFi run supplied, the audit checks up to 500 deterministically sampled
positions per group: shared discordant calls, shared agreeing heterozygotes,
asm10-only heterozygotes and asm20-only heterozygotes. Exact per-base molecule
counts and quality/observability criteria are reported. This is direct read
support, not independent sequencing validation; it uses the same reads involved
in assembly and a single-haplotype reference. BAQ is disabled for these direct
counts, explicitly distinct from production variant calling. Group sampling is
stratified, so pooled support fractions are not genome-wide error-rate estimates.
No cuts, input masks, calls or PSMC parameters are changed by the audit.

Return `comparisons/psmc-reliability-06.review.tar.gz`. It contains the Markdown
report, region/window tables, sampled molecule evidence, provenance and R figure.
Without `--hifi-run`, genotype agreement is assessed and read support is marked
unavailable. Neither agreement nor increased callable coverage alone establishes
which preset is biologically correct.

## Full 20 kb and 10 kb length trials with asm20

After syncing, from the project root:

```bash
module load miniconda3
bash gcl_genome_assembly/scripts/comparisons/submit_psmc_length_tests.sh \
  comparisons/psmc-full-presets-05 \
  comparisons/psmc-length-tests-07 \
  comparisons/psmc-hifi-comparison-03
```

The existing completed `asm20` run is the retained 50 kb control. Two fresh
Nextflow sessions run the production PSMC workflow with `asm20` and 20 kb / 10 kb
minimum assembly-alignment lengths, respectively. All other input filters,
scope, mutation/generation scaling, fitting settings, bootstrap count and seed
are cloned from that control. No assemblies or HiFi mappings are rebuilt. PAF
and SAM alignments/calls are regenerated, with matching explicit thresholds in
both mask and variant-filter recipes. The installed tools remain unchanged.
The parameter is general purpose; the chosen individual comes from source run
metadata, not a sample-specific rule.

Each trial automatically triggers a reliability audit against the 50 kb control,
using the completed HiFi BAM from the third argument. Its variant-calling step
can still be running or failed: only the successful mapping and verified matching
reference are required. Both trials and their subsequent audits can overlap.
The same seeded site sampling, genotype concordance and direct read-support
criteria are used. Reliability TSVs now use `left` / `right` columns and record
their actual run labels explicitly in `summary.json` and the Markdown report;
there is no legacy-output parser.

An after-any collector produces:

```text
comparisons/psmc-length-tests-07.review.tar.gz
```

The archive's `review/README.md` compares final callable sequence,
heterozygosity, retained scaffolds and observed PSMC histories. It links both
reliability reports. Compact provenance, generated-recipe diffs, all fit files,
input QC, per-window evidence, sampled molecule counts and R figures are retained;
large consensus, VCF, BAM and work directories stay on Crest. Failed or incomplete
workflows/audits produce an explicitly incomplete combined review rather than a
completion PASS. Trial completion does not establish biological acceptance.

The collector/controller log paths and submitted job IDs are printed by the
submission script. Per-trial outputs also remain in `asm20_20kb/` and
`asm20_10kb/` for normal Nextflow inspection. Fresh production runs still default
to asm5 and 50 kb until the experimental settings are accepted.

Native checks include a planted 30 kb / 15 kb alignment fixture: each region must
gain both confident mask coverage and its matching heterozygous call only at the
appropriate threshold. Full native PSMC workflow runs exercise both new settings.
`tests/run_psmc_length_batch_checks.py` verifies preparation, actual retained
outputs, automatic-equivalent reliability audits, R rendering, compact collection
and explicit failed-audit status.
