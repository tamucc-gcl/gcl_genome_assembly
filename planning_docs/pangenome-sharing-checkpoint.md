# CLIP sharing checkpoint

This batch adds Panacus 0.5.2 coverage histograms and expected growth tables, separately grouped by assembly haplotype and biological individual. The adapter reads path identities only; Panacus counts graph sequence. The reference and its partner retain Cactus names and are grouped via pangenome_identity.tsv. No Cactus process or input was changed.

Unknown graph paths and missing ledger members fail explicitly. Each histogram must contain coverage bins 0 through the ledger denominator, including zero-valued bins. Thus absence of universally shared sequence cannot shrink the cohort denominator. Haplotype and individual histograms must agree on graph length. No reference individual is counted twice, and no independence of paired haplotypes is assumed.

Outputs are in pangenome/<taxid>/sharing/: two histograms, two growth tables, explicit path grouping files, denominator JSON, sharing_summary.tsv, versions.tsv, and a collapsed-section Markdown report. Core/softcore/shell/cloud use the existing parameterized minima with precedence in that order. Private is separately reported because it can overlap a tier in small cohorts. Zero-coverage graph sequence is separately reported. Percentages describe distinct represented graph bp, not sequence length or repeat copy number in each assembly.

Panacus reads compressed CLIP directly. Growth uses its small histograms, avoiding another whole-graph pass. The two graph counts run sequentially within the existing 64 GB / 8 CPU / 4h Panacus allocation; adequacy must be checked on the real graph. No large node-by-path table is generated. The pangenome_growth toggle gates this batch; bp is fixed for sharing regardless of the legacy growth_count setting.

Per-chromosome and per-haplotype attribution, chromosome-by-haplotype categories, and PCA/tree remain subsequent checkpoints. Do not interpret these global histograms as those outputs. Validate tool-supported regional attribution semantics before implementation, including shared nodes across chromosomes and unplaced sequence.

## Checks and resumed run

Keep work, results and the single-library input baseline. After syncing:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv validate
```

After validation succeeds:

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv graph_build
```

Cactus should resume from cache with unchanged inputs/configuration; the sharing process is new. If Cactus unexpectedly starts reconstruction, stop and inspect the cache difference before proceeding. The new unit tests and Nextflow syntax have not been executed locally.

## Collect the checkpoint

From the cluster project root, after success:

```bash
mkdir -p comparisons
stamp=$(date +%Y%m%d-%H%M%S)
find genome_assembly/pangenome -type f \( \
  -path '*/sharing/*' -o -name 'pangenome_identity.tsv' \
  -o -name '*.vg_stats.txt' -o -name 'pangenome_artifact_checks.tsv' \
  \) -print0 > "comparisons/pangenome-sharing-${stamp}.files"
tar --null -czf "comparisons/pangenome-sharing-${stamp}.tar.gz" \
  -T "comparisons/pangenome-sharing-${stamp}.files"
ls -ltrh comparisons/pangenome-sharing-*.tar.gz logs/nextflow_graph_build_*.log
```

Send the archive plus the corresponding Nextflow and Slurm logs. Check 10 haplotypes/5 individuals for this current cohort, agreement with CLIP graph length, reference grouping, memory/time, and successful growth output. These expected counts are review criteria for this dataset, never pipeline constants.

Interface verified against toolmaker source:
- https://github.com/codialab/panacus/blob/v0.5.2/src/commands/hist.rs
- https://github.com/codialab/panacus/blob/v0.5.2/src/commands/growth.rs
- https://github.com/codialab/panacus/blob/v0.5.2/src/io.rs
