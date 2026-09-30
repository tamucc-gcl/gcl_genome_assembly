# Reference and pangenome input audit checkpoint

## Preflight repair after job 1501420

The configured image failed before parsing help because Toil could not create `/home/jselwyn/.toil` inside Singularity. The preflight label now uses Singularity's `--home` option to mount the persistent task directory as a writable container home. This change is scoped to the preflight task; it does not change assembly or graph computation. Error output remains preserved. Retry `graph_validate` after syncing; no results cleanup is needed. The fix has been statically checked but awaits cluster execution.

30 September 2026. This begins the next refactor batch, not a declaration that all remaining assembly paths or pangenome analyses are validated.

## Changes

Eligibility schema 2 records the harmonization reference, the chosen graph reference and the selection reason. A valid but no-longer-eligible reference can be replaced using the existing ranking, now explicitly reported. Duplicate, unknown and cross-species reference assignments fail rather than disappearing into fallback selection. Singleton cohorts remain withheld from graph construction while assembly reporting continues. Fragmented phased members remain eligible without becoming reference candidates.

Before Cactus, PANGENOME_INPUT_AUDIT validates biological individual count, unique assembly/graph identities, taxid agreement, reference identity and actual reference chromosome presence. It reads the input FASTAs, rejects empty/duplicate sequence IDs and records SHA256, sequence counts and lengths. Its identity ledger is published next to graph outputs and linked by the build report. This is an input integrity check, not an assembly quality certification.

The existing CLIP/GREF(CLIP)/FULL roles and Cactus 3.2.1 image selection are retained. No graph algorithm upgrade, biological threshold change, or new annotation analysis is included. Remaining graph artifact and downstream method validation is subsequent work.

## Bounded remote validation

After syncing, from the project root:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv validate
```

After that passes:

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv graph_validate
```

New mode graph_validate resumes the assembly workflow, sets run_pangenome=true and pangenome_validate_only=true, and performs eligibility, the input audit, and configured-image help/capability checks. It does not instantiate Cactus construction or pangenome analysis. It is not a standalone graph entry point: upstream tasks can run if their cache is unavailable. Preserve the work directory and prior cache. The eligibility/report tasks should update; the computational assembly stages have not been edited in this batch.

Nine new synthetic tests cover reference mismatch/replacement and graph input integrity. They have not been run locally, nor has Nextflow compilation or biological analysis. Static diff checks passed. Cluster runtime validation is required.

## Return files

After completion, package the audit outputs (the taxid below is specific to this test run, not a pipeline rule):

```bash
mkdir -p comparisons
stamp=$(date +%Y%m%d-%H%M%S)
tar -czf "comparisons/pangenome-input-audit-${stamp}.tar.gz" \
  genome_assembly/assembly_eligibility.json \
  genome_assembly/assembly_eligibility.tsv \
  genome_assembly/assembly_eligibility.md \
  genome_assembly/pangenome/373251/pangenome_identity.tsv \
  genome_assembly/pangenome/373251/pangenome_input_audit.json \
  genome_assembly/pangenome/373251/pangenome_input_report.md
ls -ltrh comparisons/pangenome-input-audit-*.tar.gz logs/nextflow_graph_validate_*.log
```

Send the new archive and corresponding Nextflow and Slurm logs. If a cohort is withheld, graph audit files will not exist; send eligibility and logs instead. Previously published graphs are not evidence of a new build in validation-only mode.

After acceptance, finish construction-output contracts before launching a long graph build, then implement Panacus-centered summaries and retained downstream analyses. Short-read and mixed-input integration validation remain outstanding. CTlk chr12 investigation remains deferred until broader improvements and the two-library rerun.
