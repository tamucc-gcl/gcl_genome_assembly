# Haplotype comparison checkpoint

## Implemented

Panacus 0.5.2 computes bp-weighted Jaccard similarity using the audited haplotype grouping from the accepted CLIP sharing task. A separate tidyverse/ggplot2/patchwork reporting task validates matrix identities, dimensions, symmetry, diagonal and finite [0,1] values, then produces a heatmap, classical PCoA, and descriptive neighbour-joining topology. It writes the distance matrix, PCoA coordinates and all eigenvalues, plus Newick/edge lengths when a tree is meaningful. No individual population analysis is added.

PCoA uses 1 minus similarity without silently correcting distances. Axis percentages use positive inertia; the negative inertia fraction is reported. A zero or rank-one matrix produces explicit absent-axis labels. Fewer than three haplotypes or all-zero distances skip NJ with an explanation. NJ topology is shown without branch-length scaling or a biological root; exact lengths including negatives remain in Newick. This avoids presenting a graph-sharing clustering as a phylogenetic result.

The existing pangenome_popstruct parameter gates this stage, within the sharing branch (pangenome_growth must also be enabled). Outputs publish to pangenome/<taxid>/comparison/. No graph construction, sharing, or accepted sharing-figure process was edited. The matrix computation uses the existing Panacus allocation, 8 CPUs/64 GB/4h; its adequacy must be assessed on the full graph. Resource usage is not yet measured. Plotting has its own pinned R environment and cache.

## Regional attribution implementation (same checkpoint)

The returned fixture proves that Panacus subset histograms do NOT preserve whole-cohort coverage: selecting A1 chr1 moves its 6 bp to coverage one, instead of retaining 2 private bp and 4 core bp. The membership export supplies node IDs and 0/1 group membership, not bp lengths, despite its bp header. Node lengths and path membership must therefore be joined separately, with duplicate traversals and walk fragments handled explicitly.

Source inspection of Panacus 0.5.2 src/analyses/table.rs also confirms that the complete text table is assembled in memory before output. A shell pipe is not a streaming-memory remedy. Do not activate a 145-million-row full export without a bounded strategy and resource check.

The new PANGENOME_NODE_COVERAGE uses Panacus table --total, exporting one global coverage count per node instead of a column per haplotype. Its compressed intermediate remains in work, not published results. Panacus still builds this smaller table in memory; this is a bounded reduction, not a claim of fully streaming execution. The existing 64 GB allocation must be checked on the real graph.

PANGENOME_REGIONAL joins these published-tool counts to GFA segment lengths and walks with NumPy arrays. It supports sparse integer segment IDs, repeated nodes, multiple W fragments of one path and P paths. It does not reconstruct haplotype coverage itself. Sorting/index arrays and temporary spooled walk text require local scratch; the adapter requests 24 GB and 24h, with actual usage to be reviewed. The routine comparison output must not include the large node table or temporary walks.

It publishes regional_coverage.tsv (exact coverage bins), regional_sharing.tsv (bp and percentages for configured tiers plus private), regional_audit.json and regional_report.md. Scopes are haplotype, chromosome, and chromosome_haplotype. Harmonized chrN_part names collapse to chrN; composite scaffolds retain separate labels; other names go to unplaced. Missing chromosome/haplotype combinations are absent rows, not biological zeros. Chromosome unions may share graph nodes and are not additive. Whole-haplotype unions deduplicate across chromosomes.

Required checks: global node coverage must exactly reconstruct the accepted Panacus histogram; summed per-haplotype bp at coverage k must equal k times that histogram bin. Tiny tests cover sparse IDs, repeat traversals, fragmented walks, cross-chromosome sharing, composite names and incorrect global coverage. Tests run inside the adapter's NumPy environment before real aggregation; local unit discovery skips those tests if NumPy is unavailable. No test or analysis was executed on the development computer. Regional plots are deliberately deferred until the returned attribution tables pass this first real-data audit.

The similarity fixture, in contrast, matches hand calculations: A1/A2=7/9, A1/B1=4/9, A2/B1=4/7. That independently supports proceeding with the small exported similarity matrix first. Source verified at https://github.com/codialab/panacus/blob/v0.5.2/src/analyses/similarity.rs and https://github.com/codialab/panacus/blob/v0.5.2/src/analyses/table.rs .

## Run and return

Static review only on the development computer; R and Nextflow were not executed here. After syncing:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv validate
# After validation succeeds:
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv graph_build
```

Expect new PANGENOME_SIMILARITY, PANGENOME_COMPARISON_PLOTS, PANGENOME_NODE_COVERAGE and PANGENOME_REGIONAL tasks plus the routine assembly summary. Cactus, sharing and sharing plots should remain cached. Keep all current work and results. If a task exceeds memory, return its logs and usage before changing settings; do not rebuild Cactus.

After success:

```bash
mkdir -p comparisons
stamp=$(date +%Y%m%d-%H%M%S)
find genome_assembly/pangenome -type f \( -path '*/comparison/*' -o -path '*/regional/*' \
  -o -name pangenome_identity.tsv \) -print0 > "comparisons/pangenome-comparison-${stamp}.files"
tar --null -czf "comparisons/pangenome-comparison-${stamp}.tar.gz" \
  -T "comparisons/pangenome-comparison-${stamp}.files"
ls -ltrh comparisons/pangenome-comparison-*.tar.gz logs/nextflow_graph_build_*.log
```

Return the archive and corresponding Nextflow/Slurm logs. To include memory measurements, use the similarity job ID printed in the Nextflow log:

```bash
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
```

Acceptance: all ten current haplotype labels present, no extra reference individual, readable heatmap/PCoA/tree, audited eigenvalues, and expected cache reuse. Counts are dataset-specific review expectations, never pipeline constants.
