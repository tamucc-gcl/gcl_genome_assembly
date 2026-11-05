# Component cleanup and retained prototype inventory — core freeze 2026-11-05

This inventory reviews 79 sources originally proposed for retirement under
`modules/` (39), `py_scripts/` (27), `r_scripts/` (9) and `workflows/` (4), including
four already-deprecated backups. At the user's request, 27 optional-analysis
sources are retained unchanged for reuse: variant classification/inversion rescue/
rearrangements, private/accessory evidence, progressive graphs and 2D layouts.
The other 52 sources remain retired. It distinguishes a working replacement from
a deferred analysis. An uncalled implementation is not evidence that its scientific purpose
was rejected, nor that a future implementation is already complete.

Original sources remain in the
[pre-cleanup revision](https://github.com/tamucc-gcl/gcl_genome_assembly/tree/341193c18878886ec18e26c3c843a44507e6162a).
The [13-category checklist](../planning_docs/remaining-pipeline-work-2026-10-15.md)
records future scope. No new analysis is enabled by this document. Retained prototypes are not included
by the production workflow; removed resource/default settings are not restored.
Their old process definitions, version records, environment requirements and
biological assumptions need review before reuse. Their presence is not a promise
that the previous flags can run them.

The restored pipeline sketches in `assets/` are retained. Updating them to the
current architecture is a separate later task; they are historical sketches until
then, not the current interface specification.

## Chimera detection and evidence

**Original use:** Enumerated AGP joins, detected chromosome-block transitions, checked coordinate provenance, lifted Hi-C pairs and gathered local read, sequence, repeat and graph evidence.

**Replaced in the core.** The maintained `workflows/chimera.nf` and `modules/misassembly.nf` use `misassembly_discover.py`, `misassembly_core.py`, `misassembly_assess.py`, `misassembly_map.py`, `misassembly_junction.py` and `misassembly_enrich.py`. Shared chimera controls, tracks, interval and graph helpers remain. Manual cutting remains supported; automatic cutting remains deferred.

Retired files:

- `py_scripts/agp_joins.py`
- `py_scripts/chimera_adjudicate.py`
- `py_scripts/chimera_blocks.py`
- `py_scripts/chimera_coordinate_guard.py`
- `py_scripts/chimera_evidence.py`
- `py_scripts/chimera_hic_pairs.py`
- `py_scripts/chimera_sequence_context.py`

## Chimera review and report layers

**Original use:** Combined earlier candidate decisions, transition families, per-assembly reports and cohort review indexes, including older HTML/SVG/Markdown presentation.

**Replaced in the core.** `misassembly_report.py`, `plot_misassembly_evidence.R`, `plot_junction_placements.R` and the main report now present the evidence and manual-cut interface. Combined evidence figures, placement figures and separate contact maps are retained.

Retired files:

- `py_scripts/chimera_decisions.py`
- `py_scripts/chimera_markdown.py`
- `py_scripts/chimera_report.py`
- `py_scripts/chimera_review.py`
- `py_scripts/chimera_review_index.py`

## Focused junction experiments

**Original use:** Focused experimental helpers for assessing selected junctions and building detailed review evidence.

**General evidence promoted into the core.** Maintained junction/enrichment helpers and R placement figures implement the accepted production analysis. Temporary investigation orchestration is not a second supported pipeline.

Retired files:

- `py_scripts/junction_assessment.py`
- `py_scripts/junction_focus.py`
- `py_scripts/junction_review.py`

## Historical join recovery and origin

**Original use:** Traced joins through earlier assembly stages and recovered older joins with unique sequence flanks for review.

**Partially superseded; historical comparison retired.** Current native-contig/graph context is retained. Recovery from older runs is not a current core requirement or a scheduled extension. This is not a claim that all historical tracing has an identical replacement.

Retired files:

- `py_scripts/chimera_origin.py`
- `py_scripts/recover_older_joins.py`

## Organelle wrappers

**Original use:** MitoHiFi-bait filtering of nuclear contigs and an older selector between HiFi and shotgun organelle assembly routes.

**Replaced in the core.** `modules/filter_organelle.nf` handles available organelle baits; `workflows/organelle.nf`, `shortread_organelle.nf` and `organelle_annotation.nf` implement the current input-aware routes and failure/status contracts.

Retired files:

- `modules/filter_mito_contigs.nf`
- `workflows/organelle_assembly.nf`

## Duplicate read-QC composition

**Original use:** Duplicate FastQC/MultiQC workflow definitions and a pairs-based Hi-C statistics/contact-map/coverage wrapper.

**Replaced in the core.** Maintained `workflows/hic_qc.nf`, `hifi_qc.nf`, `hic_mapping_qc.nf`, `hic_scaffold_qc.nf` and final Hi-C mapping/metrics modules provide current stage/final QC. Removing the old wrapper does not remove pairs-based metrics or contact maps.

Retired files:

- `modules/hic_qc.nf`
- `modules/hifi_qc.nf`
- `workflows/hic_qc_from_pairs.nf`

## Old FCS composition

**Original use:** An older combined FCS screening/cleaning and optional Blob evidence workflow.

**Split disposition.** FCS screening/cleaning remains in the maintained assembly refinement/decontamination route. The optional Blob evidence branch is deferred to category 9, not replaced by FCS.

Retired files:

- `workflows/decontamination.nf`

## Blob taxonomic/coverage evidence

**Original use:** Built and plotted Blob directories using taxonomic hits, DIAMOND databases and read coverage, with an FCS/Blob evidence summary.

**Deferred, category 9.** No complete validated replacement currently. Ordinary final read mapping is not equivalent to Blob evidence. Active database setup/build helpers and FCS remain. Rebuild evidence first; assess explicit cleaning separately.

Retired files:

- `modules/blobtools2_create.nf`
- `modules/blobtools2_viewplot.nf`
- `modules/diamond_blastx.nf`
- `modules/diamond_db_get.nf`
- `modules/fcs_blob_evidence_report.nf`
- `modules/map_reads_minimap2.nf`
- `workflows/generate_decontam_evidence.nf`

## Graph sharing and growth

**Original use:** Panacus growth/openness summaries, graph-node coverage attributed to haplotypes and private runs.

**Core summary replaced; private validation deferred.** `pangenome_sharing` and `pangenome_regional` modules/helpers and their R figures provide accepted sharing and attribution. Validated private sequence is a distinct category-11 requirement. Growth fitting and explicit open/closed summaries are now agreed follow-up work for both haplotype and individual grouping.

Retired files:

- `modules/pangenome_growth.nf`
- `modules/pangenome_hap_coverage.nf`
- `py_scripts/gfa_hap_coverage.py`

## Combined legacy graph plots

**Original use:** Rendered coverage/growth, fitted Heaps-law summaries with uncertainty, and variant class/size spectra in one plotting layer.

**Partial replacement.** Current sharing/regional R figures replace the core presentation. Interpreted SV spectra await category 12. Growth fitting/open-closed summaries are now requested for category 13. The old uncertainty formula will not be reused unchanged; it needs a covariance-aware, explicitly labelled replacement.

The historical methods and new decisions are specifically:

- **Power-law fit:** the old renderer fitted `pan(m) = k * m^gamma` by log-log
  least squares, starting at three haplotypes. Reimplementation is agreed for both haplotypes and individuals, with minimum-size and fit guards; initial delivery need not extrapolate.
- **Open/closed label:** it classified `gamma > 0.1` as open and everything else
  as closed. The follow-up restores this as an explicit configurable operational cutoff, with the exponent and sample size reported; it is not a universal saturation criterion.
- **Uncertainty:** it plotted `mean +/- 1.96 * sqrt(sum(h_i * p_i * (1-p_i)))`.
  This expression has no between-position covariance terms and does not model
  the paired sampling of haplotypes within individuals. It is not automatically
  a validated confidence interval for the biological growth curve. Choose and
  validate uncertainty estimation separately from the plotting renderer.

Current growth curves still come directly from Panacus for both haplotype and
individual grouping. The [agreed follow-up](../planning_docs/pangenome-qc-growth-followup.md) defines the added fits, classification, guards and uncertainty distinction.

Retired files:

- `modules/pangenome_plots.nf`
- `r_scripts/pangenome_plots.R`

## Graph similarity and ordination

**Original use:** Used graph-based haplotype similarity to generate ordination and neighbour-joining summaries; names sometimes said PCA where the analysis was PCoA.

**Replaced with a changed method.** `modules/pangenome_similarity.nf` and `r_scripts/pangenome_comparison_plots.R` use bp-weighted distinct-sequence Jaccard summaries, PCoA and NJ. Same descriptive purpose, not numerical equivalence or validated population-genetic inference.

Retired files:

- `modules/pangenome_pca_nj.nf`
- `modules/pangenome_popstruct.nf`
- `r_scripts/pangenome_popstruct.R`

## Input representation coverage

**Original use:** Joined Cactus sample statistics with chromosome/composite names to expose missing or weakly represented inputs.

**Purpose replaced and expanded.** Source-coordinate, sequence and clipping/build audits plus regional context now distinguish missing representation from graph/path length. These are clearer contracts, not the exact old coverage table.

Retired files:

- `modules/pangenome_input_coverage.nf`
- `py_scripts/pangenome_input_coverage.py`

## Additional graph-QC measurements and MultiQC

**Original use:** Collected acyclicity, degree/depth, tangling/linearity and alignment edit-rate summaries, and made graph/variant MultiQC sections.

**Only partially replaced.** Current core graph/source/sequence/clipping audits validate accepted construction. The missing intrinsic measurements and pangenome MultiQC HTML/data report are now agreed implementation work for categories 10/13. MultiQC will supplement and be linked from Markdown; it is no longer excluded by the reporting direction.

The correspondence is specific rather than a one-for-one replacement:

| Old measurement or output | Current state |
| --- | --- |
| Basic graph size/count statistics | Retained separately in `modules/pangenome_stats.nf`: `vg stats -lz` and `odgi stats -S` outputs. These were not removed. |
| Source/path representation and FULL/CLIP accounting | Current build, source and sequence audits provide explicit accounting, identities, source-coordinate coverage and clipping dispositions. These address construction correctness; they are not topology metrics. |
| `vg stats -A` acyclicity | Not reproduced as that explicit QC observation. Cyclicity was informational even in the old module, not an automatic defect. |
| Mean and maximum node degree | Old `odgi degree -S` parsed summaries are not reproduced as equivalent current QC fields. |
| Node-depth diagnostic | Old code retained only the first 20 lines of `odgi depth` in raw diagnostics, not a full validated genome-wide depth summary. That diagnostic is not reproduced. |
| Mean link length and sum of path-node distances | Old optimized-copy `odgi stats -l -s` measurements are not reproduced as equivalent fields. |
| Cactus GAF edit rate, identity, aligned bp and alignment count | Old pooled match/block-length calculations are not reproduced as equivalent current QC metrics. Source coverage audits measure something different. Future read-to-graph self-mapping is also a distinct analysis. |
| Whole-genome/per-chromosome MultiQC YAML | Old `odgi stats -m -s -g -d -l` export is not recreated. Its extra component acyclicity, self-loop/composition and layout summaries are not guaranteed by current basic statistics. |
| `bcftools stats` | Still produced by `modules/pangenome_variant_audit.nf` for the current Cactus-exported VCF, with sample/coordinate checks. This is not the old processed/classified catalogue or necessarily the same variant population. |
| Combined pangenome MultiQC HTML/data directory | Not currently recreated; now requested as a supplement linked from Markdown in the next pangenome work. |

Retired files:

- `modules/pangenome_qc.nf`
- `modules/multiqc_pangenome.nf`
- `modules/pangenome_odgi_stats_mqc.nf`

## Graph reuse manifest and graph-derived reference

**Original use:** Exported a graph-product reuse manifest and extracted a graph-path-named FASTA/index for graph mapping or variant workflows.

**Partial replacement/deferred exports.** Typed graph roles, provenance and published-path contracts cover current core use. Current shared-reference selection chooses a final linear assembly; it is not an equivalent graph-path FASTA export. Reassess external-consumer manifests and graph reference exports in categories 10/13.

Retired files:

- `modules/pangenome_manifest.nf`
- `modules/pangenome_ref_fasta.nf`

## Legacy pangenome report

**Original use:** Assembled older pangenome report fragments and data across graph analyses.

**Replaced in the core.** The maintained main Markdown renderer, artifact registry and linked analysis reports cover enabled core analyses. Future extensions must add their own evidence and report integration.

Retired files:

- `modules/pangenome_report.nf`
- `r_scripts/pangenome_report.R`

## Variant decomposition and classification

**Original use:** Decomposed/bubbled graph variants, split alleles and classified graph traversals by topology, frequency and size.

**Deferred, category 12.** Current Cactus exports and variant audits/counts do not replace an interpreted or finely decomposed SV catalogue. Method validation precedes reimplementation.

Retained prototype files (restored unchanged; not enabled in production):

- `modules/pangenome_variants.nf`
- `modules/pangenome_classify.nf`
- `py_scripts/classify_variants.py`

## Inversions and path rearrangements

**Original use:** Rescued potential hidden inversions by allele self-alignment; indexed/projected graph paths with odgi untangle; inferred and plotted rearrangements; supplied pilot SV diagnostics. The explicit step-index helper addressed an old odgi/Cactus indexing failure.

**Deferred, category 12.** Retained assembly-SV/SyRI code is a pilot, not a replacement for validated graph-wide rearrangement calls. Reevaluate the old version-specific indexing workaround only if the future route still needs it.

Retained prototype files (restored unchanged; not enabled in production):

- `modules/pangenome_inversion_rescue.nf`
- `py_scripts/rescue_inversions.py`
- `modules/pangenome_stepindex.nf`
- `modules/pangenome_untangle.nf`
- `modules/pangenome_rearrange.nf`
- `py_scripts/rearrange_from_untangle.py`
- `modules/pangenome_rearrange_plots.nf`
- `r_scripts/pangenome_rearrange_plots.R`
- `r_scripts/sv_pilot_diagnostics.R`

## Private/accessory validation

**Original use:** Extracted putatively private graph segments, searched competing assembly placements, measured read/k-mer support and combined evidence with controls and plots.

**Deferred, category 11.** No complete validated replacement. A cloud/private sharing count is not equivalent to demonstrated individual- or haplotype-specific biological sequence.

Retained prototype files (restored unchanged; not enabled in production):

- `modules/pangenome_private_fasta.nf`
- `modules/pangenome_private_index.nf`
- `modules/pangenome_private_map.nf`
- `modules/pangenome_private_kmer.nf`
- `modules/pangenome_private_join.nf`
- `modules/pangenome_private_plots.nf`
- `py_scripts/extract_private_fasta.py`
- `py_scripts/summarise_private_map.py`
- `py_scripts/summarise_private_kmer.py`
- `py_scripts/join_private_evidence.py`
- `r_scripts/pangenome_private_plots.R`

## Progressive graph experiments

**Original use:** Incrementally added assemblies to minigraph to investigate growth and order/reference sensitivity.

**Deferred, category 13.** Current Panacus growth is not the same experiment. Revisit utility and design with additional graph products.

Retained prototype files (restored unchanged; not enabled in production):

- `modules/pangenome_progressive.nf`
- `modules/pangenome_progressive_plot.nf`
- `r_scripts/pangenome_progressive.R`

## Two-dimensional graph layouts

**Original use:** Sorted, laid out and drew per-chromosome graph topology with odgi.

**Deferred, category 13.** Current chromosome painting, ribbons and dotplots do not replace a graph topology layout.

Retained prototype files (restored unchanged; not enabled in production):

- `modules/pangenome_2d_viz.nf`

## Report layout mockup

**Original use:** Generated synthetic presentation/layout previews during report design.

**Development-only prototype retired.** The real-data report renderer implements the accepted layout; no scientific evidence calculation was removed.

Retired files:

- `r_scripts/report_layout_preview.R`

## Already-deprecated backups

**Original use:** Old custom telomere scanner and duplicate summary-report snapshots.

**Superseded backups.** Current `modules/tidk.nf` provides telomere analysis; the maintained summary-report module/R renderer remain. No second executable backup tree is needed.

Retired files:

- `modules/deprecated/scan_telomeres.nf`
- `py_scripts/deprecated/scan_telomeres.py`
- `modules/deprecated/summary_report.nf.bak`
- `r_scripts/deprecated/generate_summary_report.R.bak`

