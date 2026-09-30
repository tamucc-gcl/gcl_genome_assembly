# First graph construction checkpoint

The input archive dated 20260930-135059 passed its input audit: ten haplotypes, five individuals, reference CMat_203_hap1 retained from harmonization, CPla included as nonvoting graph members. That archive alone did not verify Cactus preflight completion; the full build continues to require a successful capability check.

Required exports are checked in the separate graph contract task. Exact CLIP output channels are optional on the expensive Cactus task so a missing export is diagnosed after construction, rather than failing output collection there. The audit records every expected role, missing/empty status, file size and a readable standard GREF VCF header. A header-only VCF is valid; lack of variant records is not itself failure. Binary integrity, index compatibility and biological correctness still require review of the actual build. No fallback graph is substituted.

Graph roles remain CLIP biological, GREF(CLIP) variants, FULL GFA and clipping statistics QC support. Raw GREF VCF and segment map remain optional until the pinned build's exports are inspected. No annotation, new population analysis or graph algorithm upgrade is included.

## Run after syncing

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv validate
```

After successful checks, the new graph_build mode enables the real, potentially long graph run. It leaves chimera cuts and other post-assembly analyses disabled:

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv graph_build
```

Keep the same single-library baseline for this construction check. Do not delete work or results. Upstream assemblies should resume; the graph build and its checks are new work. Input/preflight tasks may reuse their successful cache. Preserve the Nextflow and Slurm logs.

## Return reports

After success, package small reports only:

```bash
mkdir -p comparisons
stamp=$(date +%Y%m%d-%H%M%S)
tar -czf "comparisons/pangenome-build-${stamp}.tar.gz" \
  genome_assembly/assembly_eligibility.json \
  genome_assembly/pangenome/373251/pangenome_identity.tsv \
  genome_assembly/pangenome/373251/pangenome_input_audit.json \
  genome_assembly/pangenome/373251/pangenome_manifest.tsv \
  genome_assembly/pangenome/373251/pangenome_artifact_checks.tsv \
  genome_assembly/pangenome/373251/pangenome_build_report.md
ls -ltrh comparisons/pangenome-build-*.tar.gz logs/nextflow_graph_build_*.log
```

Send the archive and corresponding logs. If the output audit fails, its TSV remains in its persistent task directory; send that TSV, roles.json and the error log. Do not rerun Cactus from scratch to fix a contract filename assumption. Static diff checks passed; tests and pipeline execution were not run locally.

Panacus summaries, ordination and revised structural-variant analyses remain the next implementation batch; this first run supplies construction artifacts and basic CLIP statistics only.
