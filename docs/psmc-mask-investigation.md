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
