# Full-pipeline chimera evidence validation

Use the normal pipeline entrypoint, not the isolated manifest harness. The harness assumed a retained reference PAF for each assembly and failed when the reference assembly had no such published file. The full workflow produces its own inputs and supplies HiFi reads through the normal sample channels.

After syncing the updated repository, submit from the Crest project directory:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch
```

The launcher delegates to run_assembly_checkpoint.sbatch in assembly mode: main.nf, -profile slurm, -resume, the main launch/cache directory `.nf/assembly`, and the main `work` directory. Existing checkpoint parameters are preserved. Worker processes use the standard per-label CPU, memory, time and conda environment definitions; the launcher itself requests one CPU and 8 GB for orchestration.

Default inputs are `data/assembly_samplesheet.csv` and `data/hic_readsets.csv`. All sample-sheet assemblies enter the normal workflow, including candidate detection, evidence, finishing and reporting. Defaults from the main checkpoint launch still govern optional stages; this launcher does not enable every optional analysis.

The third positional argument can choose a specific separate output root:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv comparisons/chimera-full-pipeline-review
```

With no third argument, outputs go to `comparisons/chimera-pipeline-JOBID/results`. Publication uses copy mode, so this run does not overwrite the main published results. A sibling `.review.tar.gz` contains the chimera evidence/report tree, main Markdown report if produced, trace, log and exit status; large sequences and BAM/SAM files are excluded.

Unchanged tasks can reuse the main cache when their source inputs, process definitions, environments and work outputs still match. Changing outdir alone should not rerun an assembler whose task does not consume that parameter. Changes to metadata, inputs, tool environments, task scripts or missing work files can invalidate cache entries. This launcher resumes the most recent session in the main assembly launch directory; do not launch another main assembly driver concurrently.

The cohort report is `results/assembly/chimeras/review/README.md`. It includes every assessed assembly, candidate ranges, readable IDs, chromosome relationships, evidence for and against cutting, observability limits, exact verified gap cuts and instructions for user-added cuts. Main-report integration of this new detailed report remains a separate update.

Local tests cover report generation and source identity checks. SLURM/Nextflow execution and cache reuse need verification from the Crest trace. Inspect CACHED entries for hifiasm; do not infer cache reuse solely from using -resume.
