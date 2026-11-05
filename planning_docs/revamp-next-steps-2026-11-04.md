# Current revamp status and next native reporting test

Latest update: the local detailed-figure batch is ready for one cached Crest test.
Use [the detailed-report acceptance notes and current command](detailed-report-local-2026-11-05.md).
The main-report-02 instructions below are retained as the prior checkpoint.

This is the current summary of the canonical 13-category checklist in
[remaining-pipeline-work-2026-10-15.md](remaining-pipeline-work-2026-10-15.md).
Historical open/closed statements in older checkpoint narratives do not supersede
the accepted category statuses at the top of that checklist.

Update 2026-11-05: main-report-02 computation and retained evidence checks passed;
the launcher failed on a valid directory link rejected by report validation. The
validator is fixed locally, and the native postchecks still need to be run. See
[the assessment and commands](main-report-02-assessment.md). Detailed figure work
has started with retained HiFi competitive-placement evidence rendered in R.

## Next Crest run: main-report-02

Sync all updated pipeline files, including the new environments and report helpers.
Run from the project directory; do not launch another controller concurrently into
the shared `.nf/assembly` session.

```bash
cd /work/birdlab/GCL/spratelloides_delicatulus_genome
mkdir -p logs

CORE_INTEGRATION_ALLOW_REBUILD=false \
sbatch gcl_genome_assembly/scripts/comparisons/run_core_integration.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/main-report-01/results \
  data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv \
  comparisons/main-report-02 \
  --run_report true --report_chromosome_painting true \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json"
```

The launcher uses normal `-resume`, the existing work directory, all-stage QC,
the accepted manual cuts, core graph and within-individual R comparisons. Publication
uses symlinks. The retained BAM manifest is the same test-reuse override recorded in
main-report-01; it is not a requirement for a fresh production run.

Assembly, scaffolding, core graph, painting alignments and compatible mappings should
remain cached. Updated report generation, chimera decision presentation, provenance
and the small annotation feature audits can execute. The integration audit requires
Cactus construction to remain cached and final assembly/reference/name-map identities
to match the accepted baseline. It does not prevent Nextflow from submitting an
unexpected task before the audit; inspect the controller's task list if reuse differs.

Post-publication validation checks expected report outputs, registered files and local
file/heading links after Nextflow exits. Results are `results/report_validation.md/json`.
`core_integration_check.md/json` records the core/cache acceptance. These checks do not
turn graph representation flags into biological failures or new cut recommendations.

Review `results/assembly_report.md`, the linked details and the five explicitly labelled
pangenome graph entries. Pla scope limitations should appear only in assessment coverage.
The rejected fragmented-assembly summary plot must not be added.
Return `comparisons/main-report-02.review.tar.gz`.

## Remaining core reporting gate: category 6

1. Evaluate this native test: successful annotation/report tasks, unchanged assembly
   and graph evidence, correct reports/links, and ordinary cache reuse.
2. Finish the detailed visualization pass. In particular, show apparent versus credible
   HiFi spans and alternate placements from the retained evidence; migrate maintained
   custom plots to R/ggplot2 where practical. Inspect Hi-C/contact, telomere and graph
   figures for readable scales, legends, captions and missing-data treatment. Existing
   native tool figures can remain where appropriate; do not claim everything is R-only.
3. Finish presentation/navigation for the existing core pangenome summaries, keeping
   technical validation accessible and the five known graph flags distinct from assembly
   defects. Richer biological presentation and post-pangenome supplements remain in the
   extension phase. The accepted, documented graph flags do not themselves reopen core
   acceptance; a concrete newly demonstrated correctness defect would.
4. Confirm report behavior for supported input/disabled/skip/failure routes, remove
   superseded active report paths/outputs after replacements pass, document and freeze
   the core parameter/output interfaces, then mark category 6 accepted.

Acceptance of this test alone does not finish every detailed-renderer or extension item.
Categories 1–5 and 7 remain accepted unless a concrete regression is found; do not
restart those investigations merely because the report presentation changed.

## Extension phase after the core gate

| Category | Remaining work |
| --- | --- |
| 8. PSMC — first extension | Both agreed diploid modes: paired same-individual haplotypes, and one assembly plus complete HiFi reads. Validate callability/masks, fitting, scheduled bootstraps, scaling and report integration, independently of graph enablement. |
| 9. BlobToolKit evidence and cleaning | Rebuild the optional evidence/report route. Validate explicit cleaning separately, including downstream invalidation if sequence genuinely changes. Existing production FCS screening is not being declared broken. |
| 10. Graph self-mapping QC | Bounded reproducible read subsets, a supported mapper/index route, standard metrics and reporting. Linear self-mapping does not substitute for graph mapping. |
| 11. Private/accessory evidence | Repair or replace the requested evidence workflow; distinguish individual/haplotype privacy, missing evidence and detection power. Decide long-term maintenance after its requested use. |
| 12. Graph SV interpretation and supplementary SyRI | Validate supported event classes, coordinates and missingness with known-answer tests and a bounded biological pilot; integrate the accepted interpretations and supplementary comparisons. |
| 13. Additional graph products and presentation | Inventory existing outputs, add only agreed useful supplements/exports/figures, and expand the biological pangenome report accordingly. |

Each extension remains independently enabled and must contribute to the report.
Independent bounded branches can run in parallel once their input contracts are ready;
do not run concurrent production resumes in the shared launch session.

Automatic chimera cutting and biological purging optimization remain deferred. They
are not requirements for finishing the current core report. New annotation/population
analysis or other previously excluded scope is not silently added to this roadmap.
