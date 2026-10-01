# Pangenome figures and regional sharing

## Accepted checkpoint and next implementation

Run 1502809 passed: 303 cached tasks, two new tasks, no failures or retries. Cactus was cached. CLIP sharing accounts for 1,949,595,189 bp, ten haplotypes and five individuals. Both growth curves cover the full cohort and their endpoints reconcile within one bp. The zero-only growth failure is resolved by supplying GFA directly.

The first figure task is now separate from computation. It plots existing Panacus output as growth curves, coverage histograms and sharing partitions, with haplotype and individual panels, PNG for GitHub and PDF for export. It publishes a collapsed Markdown fragment with relative links. It does not rerun Panacus to make plots, fit openness, or invent uncertainty bands. Keep the accepted sharing process unchanged while adding reporting.

## Inventory of original plots and proposed disposition

Source inventory: R_scripts/pangenome_plots.R, pangenome_private_plots.R, pangenome_popstruct.R, pangenome_progressive.R and pangenome_rearrange_plots.R, plus references in pangenome_report.R. This is an inventory of retained legacy source, not a claim that every figure was produced in a historical run.

| Original plot family | Decision and reason |
|---|---|
| Growth/core curves, coverage histogram, sharing partition | Retain; plot Panacus results directly. Implemented in the new figure task. |
| Heaps-law fit and open/closed interpretation | Omit from the general default; this graph/cohort does not establish such a classification. |
| Progressive construction growth | Omit duplicate default; use expected accumulation from the final CLIP graph. |
| Private bp/fraction by haplotype | Retain within a complete sharing composition figure; private-only ranking hides other categories. |
| Tier/private by chromosome | Required, extend to chromosome-by-haplotype composition with explicit denominators and unplaced sequence. |
| FULL-versus-CLIP private sequence | Compact clipping QC/sensitivity only, not a second full suite of biological plots. |
| Sharing-segment length spectra | Defer until a concrete question requires them; segments are not automatically SV events. |
| Private/control, copy ratio and special evidence panels | Exclude from standard report; retain only targeted diagnostics that have a supported use. |
| Haplotype ordination and NJ tree | Retain after validating the published-tool matrix. Label PCoA versus PCA accurately; report negative eigenvalues if present. A similarity tree is descriptive, not an inferred species history. |
| Individual ordination/tree | Excluded from this assembly pipeline, per agreed scope. |
| SV counts/size spectra, variant classes, chromosome reference footprint | Retain in the later GREF variant checkpoint, with event/allele/bp and reference-coordinate distinctions explicit. |
| Inversion, duplication, translocation and orientation panels | Retain the capability; finalize plots after the planned method comparison, not from unverified old calls. |

Useful additions: a haplotype similarity heatmap beside ordination; chromosome-by-haplotype sharing-fraction heatmaps; and a compact input-versus-retained-sequence QC view to help distinguish biological sharing patterns from clipping/fragmentation. Integrate these into the same sections, avoiding redundant figures. Absolute bp and fractions should be available together for composition plots.

## Regional accounting contract

Sharing coverage is always computed against the complete eligible cohort. Selecting one haplotype or chromosome for attribution must not redefine core as shared by only that subset. Count graph sequence once per reporting unit, distinguish this from assembled copy length, and keep unplaced paths explicit. A node found on multiple chromosomes can contribute to more than one chromosome's distinct-sequence total, so these totals must not be silently summed as a genome partition. Do not infer homology purely from similarly named unplaced/composite contigs.

The Panacus table interface supports custom grouping, but its node-by-group export can be enormous on the present 145-million-node graph. Do not enable that export on the full dataset before assessing format, resource use and a streaming/partitioned strategy. The CoverageLine source is a histogram analysis, not automatically a per-chromosome track. The tiny comparison probe tests regional subset semantics, membership export and bp-weighted similarity on hand-checkable data. This is required before selecting the adapter; it avoids another plausible-looking but incorrect tool integration.

Toolmaker interfaces inspected at tag v0.5.2:
- https://github.com/codialab/panacus/blob/v0.5.2/src/commands/table.rs
- https://github.com/codialab/panacus/blob/v0.5.2/src/commands/similarity.rs
- https://github.com/codialab/panacus/blob/v0.5.2/src/analyses/coverage_line.rs

## User-run checkpoint

No pipeline, R, Python tests or analysis was executed locally. After syncing, use the normal unit tests and validation, then the same graph_build command. The new figure task should run while sharing and Cactus remain cached. Inspect PNGs/PDFs for readable legends and correct thresholds before accepting presentation.

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv validate
# After validation succeeds:
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv graph_build
```

Run the tiny method probe using the Panacus environment recorded in the successful log (this path is specific to the current cluster, not embedded in pipeline code):

```bash
mkdir -p comparisons
sbatch --partition=normal --time=00:15:00 --mem=2G --cpus-per-task=1 \
  --output=logs/panacus-methods-%j.out \
  --wrap='module load miniconda3; conda run --no-capture-output -p /work/birdlab/.conda_builds/env-4e35a1c4ad242753bb37f95cd57242e7 bash gcl_genome_assembly/scripts/comparisons/check_panacus_regional.sh comparisons/panacus-methods-$(date +%Y%m%d-%H%M%S)'
```

After both jobs finish, collect figures and tiny probe outputs (no graph files):

```bash
stamp=$(date +%Y%m%d-%H%M%S)
find genome_assembly/pangenome -type f \( -path '*/figures/*' -o -path '*/sharing/*' \) \
  -print0 > "comparisons/pangenome-figures-${stamp}.files"
tar --null -czf "comparisons/pangenome-figures-${stamp}.tar.gz" \
  -T "comparisons/pangenome-figures-${stamp}.files"
tar -czf "comparisons/panacus-methods-${stamp}.tar.gz" comparisons/panacus-methods-*/
```

Return both archives and the corresponding run/probe logs. Keep work and results. Regional attribution and ordination are not yet wired into the production workflow; the returned tiny fixture determines the next implementation without a large speculative rerun.
