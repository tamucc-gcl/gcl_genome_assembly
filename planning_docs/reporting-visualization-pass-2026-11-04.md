# Reporting and visualization pass — 2026-11-04

Category 6 is active after core-integration-02 acceptance. Evidence calculations,
assembly decisions and graph construction are accepted inputs to this work.

## Presentation decisions from the user

- Concise main GitHub Markdown report, with key figures and linked detailed reports.
- Open with assembly quality and chromosome structure, rather than execution history.
- Show individual quality metrics: BUSCO, QV, k-mer completeness and contiguity.
- Show examples of all three proposed chromosome-structure views before choosing:
  scaffold lengths/names, chromosome painting, and a representative ribbon plot.
- Subsequent choice: A (scaffold lengths/names) and B (measured chromosome painting)
  together in the main report; C (ribbon comparisons) in details.
- Clarification: combine A+B on the same bars, not as separate figures. Each scaffold
  bar has its name, a length scale and internally painted chromosome blocks. The
  combined fictional mock-up is rendered as structure_AB_combined_concept.png and
  replaces the separate A/B figures in the main preview. Production blocks must
  come from coordinate-matched evidence; fictional boundaries are never production input.
- Quality decision: drop the summary quality figure entirely. Use a native Markdown
  table, one individual per row, with each metric cell showing haplotype 1 / haplotype 2.
  Do not render this table as an image. The preview no longer generates quality plots.
- New/revised graphics use R/ggplot2 where practical, with pinned dependencies and
  actual rendering checks. Evidence generation stays separate from presentation.

## First reviewable prototype

`r_scripts/report_layout_preview.R` reads retained accepted outputs and writes
`report-layout-preview.md` plus `reports/layout-preview/` under a supplied results
directory. It is a presentation prototype, not the production report entry point.
The original assembly_report.md is preserved.

The preview uses real final QC values, all ten final scaffold name maps, verified
cut records, reviewer dispositions and an existing R/gggenomes ribbon figure.
The painting example alone is fictional and prominently labeled in both figure
and report. Full final-coordinate painting blocks are absent from this review
archive; chromosome names must not be used to invent internal boundaries. If a
measured painting is selected, explicitly assess compatible retained alignments
and any additional mapping/liftover needed before changing the evidence inputs.

Local QA: 40 final quality points; six original transition events resolve to four
verified applied cuts and two explicit RETAIN decisions. Parent PENDING rows and
unused alternative cut locations do not create new unresolved events when a
verified child action already resolves that event. All local links in the preview
and its decision detail page resolve. Real PNGs were visually inspected; PDF font
handling and the quality legend were corrected after rendering. Tested environment:
R 4.3, tidyverse 2.0.0, ggplot2 3.5.2; session details are retained with the preview.

Pending preference: whether k-mer completeness should show haplotype values, combined
diploid values, or both with visibly distinct measurement units. The current
prototype shows haplotype values and explicitly labels that scope.

## Implementation sequence

1. Settle opening layout using rendered examples and user feedback.
2. Define a run-scoped artifact registry with published relative locations,
   coordinate stages, measurement units and availability. Consume explicit current
   workflow channels; do not scan a potentially stale publication directory in
   production. Preserve each imported report's directory when rebasing links.
3. Replace concatenation with a concise report plus linked details; fix section
   hierarchy, enabled-section omission, provenance identity and all relative links.
   Keep interactive QC ranking HTML as the agreed exception to Markdown reports.
4. Report detected events separately from human decisions and verified actions.
   Join decisions by assembly, assessed input fingerprint and event identity;
   preserve manual additions and original-to-final coordinate records. Never infer
   an applied cut merely from a checked spreadsheet row.
5. Migrate active Python evidence renderers to R using retained measurement tables,
   including alternate HiFi placements, calibrated/observable contact evidence,
   telomeres and clear cut/range annotations. Keep Hi-C contact maps separate where
   already agreed. Do not reactivate dormant renderers solely for migration.
6. Review graph sharing/representation figures and final Hi-C/QC figures for legends,
   scales, missing-data treatment, composite labels and readability. Avoid changing
   biological filtering merely to make a cleaner picture.
7. Add meaningful enabled/disabled, first-pass/reviewed/manual-addition and multi-species
   report cases; validate figure values against retained evidence and check rendered
   figures and local links. Run the real renderer locally on accepted outputs first.
8. Run a normal cached native report integration check. Accepted assembly/graph
   producers should remain cached unless a separately justified evidence change
   is explicitly chosen. No new full assembly/graph run is requested for this prototype.

PSMC remains category 8 and the first extension after the core reporting gate.

## First production implementation (native acceptance pending)

The main report now runs independently of QC/post-assembly enablement under
`run_report`. The existing detailed report is preserved as `assembly_details.md`;
its detailed visual/content redesign is deferred until the main report is accepted.
The main report uses explicit artifact channels and includes the approved Markdown
quality table, combined length/painting view, verified event outcomes, sample
outcomes/files, graph review notes/links, organelle outcomes and evidence/provenance
links. Graph fragments retained in the old detailed report have their links rebased.

New reporting evidence is explicit: `REPORT_PAINTING` maps each finalized assembly
against its species' selected final shared reference, using the existing comparison
preset/MAPQ/chain-length policy and primary mapping chains (no base-level alignment).
Self-reference uses identity labels; composite/unplaced reference labels cannot
paint an internal chromosome. The compact projection partitions query coordinates
exactly once: same-label overlaps coalesce, conflicting labels become ambiguous,
uncovered spans remain unassigned. MAPQ 255 is unavailable. This does not feed
assembly naming, chimera scoring/cutting, or graph construction.

All scaffolds remain in the evidence table. Only the longest configured count is
displayed, with displayed/total scaffold counts and sequence fractions retained.
The renderer fails on missing expected painting outputs or dictionary mismatches.
The source PAF/log/metadata and versions are retained separately from figure files.
R/ggplot2 dependencies are pinned in environments/reporting.yaml; tasks do not
install packages at runtime.

Verification so far: four painting accounting/dictionary tests; eight real-R cases
(reviewed, fresh pending, wrong fingerprint, manual additions, no organelle recovery,
primary/collapsed outputs, disabled inputs and accepted real measurements);
every reported QC measurement is independently reconciled with its source CSV.
Actual Nextflow production-workflow execution passed with
painting enabled, disabled, and multiple species using a controlled mapper and
real projector/R renderer. The full main.nf compiles and its input-only checkpoint
passes. Graph fragment rebasing/collection passes an actual Nextflow case. Native
mapping on the real final assemblies and the resulting combined figure still need
Crest review. Core assembly/graph producers should stay cached in that run.

The integration audit additionally checks main-report membership, expected painting
membership and existence of every registered report artifact. The first reporting
run should use the accepted core-integration-02 baseline in cache-only mode.

After syncing, run from the Crest project directory:

```bash
CORE_INTEGRATION_ALLOW_REBUILD=false \
sbatch gcl_genome_assembly/scripts/comparisons/run_core_integration.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-integration-02/results \
  data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv \
  comparisons/main-report-01 \
  --run_report true --report_chromosome_painting true \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json"
```

Return `comparisons/main-report-01.review.tar.gz`. New mapping is final assemblies
to the selected reference for painting (nine query mappings plus one identity
projection for this cohort), followed by report/provenance tasks. Assembly, read
mapping, chimera assessment and Cactus construction should reuse normal cache.
The audit rejects an unexpected Cactus rebuild. Existing detailed reports are
retained for the subsequent functionality/presentation review; this does not
claim migration of all active Python renderers to R.

## Subsequent report review

The native main-report-01 run has now been inspected and its execution/evidence
accounting accepted. Later presentation revisions, the selected separate-scale
unplaced summary, individual-based tables, focused detailed reports, annotation
feature-count audit and validation results are recorded in
[main-report-01-assessment.md](main-report-01-assessment.md). That assessment supersedes
the earlier "still need Crest review" status for the first report run. The latest
reporting revisions have passed local real-data and Nextflow checks and await their
next native integration run; the broader reporting category remains open for review.
