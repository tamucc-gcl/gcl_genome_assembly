# Remaining pipeline rework — reconciled 2026-10-15

This inventory reconciles the agreed implementation plan, the 2026-10-02 ordering update,
PSMC design, current active workflows and the completed Crest QC evidence. It supersedes
abbreviated remaining-work summaries, without removing agreed extensions or changing their
ordering. Optional means independently enabled per run, not forgotten from this roadmap.

## Current accepted checkpoint

- General cohort chimera evidence, explicit manual cuts and downstream finishing have been
  exercised on the real ten-haplotype cohort. Automatic cutting remains explicitly deferred.
- Chromosome-only telomere extension, record-order restoration and all-stage QC passed the
  latest uploaded real run: 116 completed tasks, 661 cached, no failures/retries. Both reviewed
  CTlk hap1 cuts remain verified. Purging/biological optimization is deferred by the user.
- Final Hi-C mapping/filtering/metrics/contact maps have now been decoupled from
  run_post_assembly and tested locally in final_only/all_stages, QC-off and plots-off cases.
  Their real Crest acceptance is pending. Local Nextflow wiring tests used tool stubs.
- Graph admission/construction, role audits and sharing/regional/similarity modules exist.
  Presence of source and synthetic tests is not acceptance of the updated real cohort.
- The always-on assembly summary and interactive QC viewer exist. The consolidated Markdown
  report and final chimera decision-status integration remain unfinished.

## Core work still required

| Workstream | Remaining work | Acceptance |
|---|---|---|
| Final Hi-C QC | Run final-coordinate mapping/filtering, retention, trans/cis and shared contact maps; confirm library handling and report inclusion | Correct final sample/haplotype counts, final-stage statistics and map artifacts, independent of post-assembly switch |
| Supported input routes | Reconcile existing tests with singleton short-read, HiFi-only, mixed read types/species, haploid/diploid representations, Hi-C-only skips and missing/duplicate inputs; execute missing bounded cases | Correct routes, identities, reasons, no cross-species/sample leakage, no silent disappearing outputs |
| Failure and organelle handling | Audit/complete required-task failure isolation, species-cohort withholding, organelle no-recovery versus tool failure, screening-reference fallback and conservative separation | Failed required members never silently shrink a graph; valid auxiliary no-result outcomes preserve nuclear assemblies; screening incompleteness is explicit |
| Shared references and provenance | Separate species linear-reference selection from graph eligibility; support singleton diploid reference selection, explicit overrides and graph-on/off consistency; record final FASTA identity and mapping compatibility | Pangenome and PSMC consume the same eligible-cohort linear reference; singleton analyses can select a reference without relaxing graph admission |
| Core pangenome | Input/image preflight; build and verify CLIP biological graph, GREF(CLIP) variant coordinates/exports, compact FULL/clipping QC; reconcile sample/path/ploidy identities | Real-cohort graph contract and variant audits pass; synthetic reference paths do not inflate biological counts |
| Required graph summaries | Verify Panacus sharing/growth and all chromosome/haplotype/chromosome-by-haplotype attributions; descriptive similarity heatmap/ordination/tree; explicit denominators and unplaced/composite sequence | Required outputs exist, known-answer tests pass and totals/overlap rules reconcile; implement any actual gaps rather than relabelling them optional |
| Consolidated report | One current-run GitHub Markdown report with essential tables/plots, collapsed sections, relative links, clear statuses; disabled sections absent; retain interactive assembly ranking as supplementary | Singleton/mixed/skip/failure/extension-disabled runs have meaningful reports; no stale directory-wide discovery; cuts and retained joins show final reviewer disposition |
| Operational cleanup | Complete tool/container/database pinning and resource/retry/disk audit; verify dependency selectors and optional channels; remove superseded code, parameters, environments and unused published outputs after replacements pass | Stable resume and bounded resource use; failed deterministic commands are not endlessly retried; user-selected output directories and cache retention preserved |

Older checkpoint documents list some of these as outstanding despite subsequent repairs.
Re-audit current code before making changes. In particular, coordinate-safe manual cuts and
chromosome-scope handling are no longer wholly unimplemented workstreams.

## Agreed extensions still to deliver after the core report/cleanup gate

| Extension | Current state | Remaining deliverable |
|---|---|---|
| PSMC | Planned; no production workflow, module or run_psmc implementation found | Independent optional PSMC workflow supporting the two agreed diploid input modes, callability, fits, scheduled bootstraps, provenance and report integration |
| BlobToolKit/decontamination evidence | Legacy modules exist; active evidence option is explicitly rejected as unsupported | Rebuild supported evidence/dataset/static-report route, then separately test explicit cleaning, removed sequence, coordinate changes and downstream invalidation |
| Graph self-mapping QC | No accepted production implementation established | Select published mapper/index workflow; reproducible bounded subsets of contributing reads, standard metrics, index provenance and report contribution |
| Temporary private/accessory evidence | Legacy modules/scripts remain; presence is not acceptance | Repair/replace identity, alignment, coverage/k-mer and control logic; distinguish private to individual/haplotype, missing evidence and detection power; decide long-term retention after requested use |
| Graph SV interpretation | Existing graph variant exports do not settle event-class interpretation | Known-answer INV/DUP/translocation/nested/repeat tests, published-method evaluation, supported classification and explicit missingness/coordinates/read support |
| Supplementary SyRI | Optional workflow/pilot exists; prior delta-path fix lacked completed acceptance in the reviewed checkpoint | Retain pilot; bounded validation and report integration if used; never substitute restricted pairwise results for graph-wide capability |
| Additional justified graph outputs | Existing products must first be inventoried | Only distinct useful additions such as focused regional diagnostics or consumer-required export/decomposition; required core summaries remain core |

Every extension must have an independent default-off switch, correct unsupported/failed
outcomes, bounded resources, stable expensive intermediates and a contribution to the existing
Markdown report. The unfinished core report is not a reason to deliver extensions without reports.

## PSMC detail that must not be dropped

1. Independent of run_pangenome, qc_mode and optional visualization toggles. A species graph
   is not a prerequisite. True biological haploids are an explicit unsupported case.
2. Two same-individual diploid haplotypes: Dipcall route with its confident-region mask.
3. One assembly representation of a diploid individual plus that individual's complete HiFi
   reads: mapping/diploid calling route. Caller/model/filter selection still needs validation;
   do not assume a human-trained model or haplotype-assigned reads are suitable.
4. Both modes use the shared final linear reference, never a synthetic GREF path. Preserve
   callable homozygous, heterozygous and unknown sites distinctly, plus user interval/mask scope.
5. Separate reference/mapping/calling, PSMC input, observed fit, bootstrap tasks and rendering.
   Plot rescaling must reuse fits. Schedule replicates through Nextflow with a concurrency cap.
6. Record callability, heterozygosity, missingness, fitting/bin parameters, input/reference
   identities and raw fits. Scaled years/Ne plots require supplied species-specific mutation
   rate and generation time; otherwise output explicitly unscaled results.
7. Accept both modes on bounded fixtures and a representative biological pilot before expanding
   to cohort/bootstrap production. Validate singleton, graph-off, mixed species, poor callability,
   truncated fits, bootstrap completeness and rescaling.

## Execution order and concurrency

1. Upcoming combined Crest run: final Hi-C QC plus validation-only pangenome preflight/input
   audit. This does not build a graph or start PSMC. The review archive includes small graph
   input audits and final Hi-C plot/statistic products.
2. While Crest runs: audit/test remaining input, failure, reference and operational contracts
   locally; prepare report interfaces. Do not sync changes to source used by live tasks.
3. After preflight acceptance: run the core graph and dependent core summaries. Independent
   graph summaries/chromosome jobs can run concurrently through the one pipeline controller.
   Develop/verify report rendering against fixtures and existing artifacts during this work.
4. Accept core graph/report/cleanup and the supported-input matrix. That finishes the core,
   not the full agreed roadmap.
5. Extension phase: PSMC first, then BlobToolKit evidence, bounded graph self-mapping,
   requested private/accessory evidence, focused graph SV/SyRI work and justified extra exports.
   Independent branches can run together once their inputs/contracts are accepted. Never run
   concurrent resumes into the shared .nf/assembly session; use one combined pipeline or
   deliberately isolated validation jobs.

Biological purging optimization remains deferred. Automatic chimera cutting is a separate
explicitly deferred revisit, not a core or extension release requirement. Gene/repeat annotation,
general population analysis, ONT/TellSeq and a standalone analysis-only entry point remain
outside this agreed rework unless separately requested.

## Source plans

- [Agreed assembly/pangenome decisions](assembly-pangenome-agreed-plan.md)
- [Authoritative post-core order and extension scope](post-core-optional-analyses-plan.md)
- [Detailed PSMC design](completed-run-review-and-psmc-plan.md)
- [Original architecture and batches](pipeline-rework-roadmap.md)
- [Optional SyRI checkpoint](assembly-sv-integration-checkpoint.md)
