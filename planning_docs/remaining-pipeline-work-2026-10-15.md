# Remaining pipeline rework — reconciled 2026-10-15

This inventory reconciles the agreed implementation plan, the 2026-10-02 ordering update,
PSMC design, current active workflows and the completed Crest QC evidence. It supersedes
abbreviated remaining-work summaries, without removing agreed extensions or changing their
ordering. Optional means independently enabled per run, not forgotten from this roadmap.

## Completion checklist — the 13 agreed categories

These IDs are the canonical progress checklist. Check a category only after its agreed
implementation and applicable acceptance evidence are complete; source presence alone is
insufficient. Update its evidence/status here when checking it off. Existing completed
substeps remain recorded below and should not be repeated without a concrete reason.

- [x] **1. Assembly/QC acceptance** — accepted 2026-10-16 from Crest run 1535116:
  55 completed tasks, 774 cached, zero failures/retries. All ten finalized assemblies have
  final Hi-C mapping/filtering/metrics and 40 contact-map plots. Final CSV/HTML retention and
  trans/cis values reconcile with pair counts and map statistics. Both reviewed cuts remain
  verified; final name maps are unchanged. Detailed evidence: local
  `comparisons/assembly-all-stage-qc-final-hic-inspection/QC-review.md` and uploaded
  `assembly-all-stage-qc.review.tar.gz`. Standalone metric-table publication was fixed/tested
  locally for the next normal resume (category 7); crowded map labels remain category 6.
- [x] **2. General-purpose input handling** — accepted 2026-10-28 from the supported
  parser/species/skip matrix, native HiFi/shotgun runs, final inventory checks and
  [native name/taxid acceptance](core-identity-reference-01-assessment.md). Name-only
  and numeric identities resolve without silently dropping samples. Final core
  integration remains category 7.
- [x] **3. Failure and organelle handling** — accepted 2026-10-28 from controlled
  positive/negative recovery/context/database/annotation tests and native HiFi,
  fragmented shotgun annotation and reference retrieval. Partial recovery remains
  explicit. Broader pipeline retry/resource/dependency cleanup remains category 7;
  biological completeness is not claimed by execution acceptance.
- [x] **4. Shared references and provenance** — reference selection is independent of graph
  enablement/admission, with singleton/override/species-isolation tests and final FASTA/FAI
  identity checks. Real run 1535231 confirms the shared reference checksum agrees with graph
  inputs. PSMC consumption belongs to category 8; graph tool-version labeling remains category 7.
- [x] **5. Core pangenome** — accepted 2026-10-25 from
  [production-02 acceptance](core-graph-production-02-assessment.md): 17 completed,
  840 cached, zero failures/retries. Both missing CTlk chr7 arms recover; CLIP
  retains 96.82%/95.88% of the composite inputs. All 6,182 sources are accounted
  through partitioning and FULL/CLIP, with no assignment overlaps or stage-coordinate
  discrepancies. Five pre-existing low-retention flags remain explicit for later
  graph QC/interpretation; they are not new regressions. Temporary diagnostics and
  report/resource cleanup remain categories 6/7. Earlier checkpoints follow.
  Construction, role checks, sharing/growth, regional attribution,
  similarity/ordination/tree and standard variant-export/count audits passed real run 1535231
  (15 completed, 826 cached; no failures/retries). Ten assemblies/five individuals; reference
  Sde-CMat_203_hap1. Finish compact FULL/clipping accounting before final category acceptance.
  Core-graph-02 completed with 4 new/839 cached tasks and no failures/retries. Its raw
  Cactus bundle exposed representation-summary parser errors (vg phase suffixes,
  synthetic editing records and BED interval interpretation), now corrected locally.
  Corrected comparison: all 4,663 measured FULL source path names survive in CLIP;
  CTlk hap1 chr4/chr14 is already absent under its expected name in FULL, requiring
  an upstream mapping/partitioning explanation before acceptance. Review:
  `comparisons/core-graph-02-inspection/QC-review.md`.
  A production `PANGENOME_SOURCE_AUDIT` now checks input FAI dictionaries against
  cached GAF/PAF and FULL/CLIP paths for every cohort member. Its Crest acceptance
  is pending: [bounded source-audit run and completion gate](core-graph-source-acceptance.md).
  Detailed original review: `comparisons/core-graph-01-inspection/QC-review.md`. No graph rebuild is
  indicated. Resource requests, Cactus/Toil version labeling and review-archive omissions
  were identified for category 7; consolidated presentation remains category 6.
  **2026-10-23 checkpoint:** v8 rebuilt successfully; the separated chr4/chr14 paths
  now survive. The retained chr7/chr12 composites retain only about half their source
  path lengths. Name presence alone had hidden this from the PASS audit. Core
  acceptance remains open pending the [cached mapping/filter/assignment coverage
  trace](core-graph-coverage-trace.md); no additional graph rebuild or assembly cuts
  are indicated before that evidence.
  The filter log subsequently confirmed whole-arm removal by the upstream deletion
  filter. A removable, pinned-tool [chromosome-isolated filtering comparison](core-graph-deletion-filter-test.md)
  is now wired into the pipeline: native, guarded and disabled filtering, each
  followed by whole-query/regional assignment. Production graph policy remains
  unchanged pending the Crest comparison and final graph acceptance.
    The deletion comparison reproduced the native exported PAF byte for byte and
    confirmed recovery under chromosome-isolated filtering. A separate upstream
    regional grouping comparison still discards an arm. The next
    [pinned regional correction test](core-graph-regional-assignment-test.md) changes
    that one comparison, audits source-coordinate ownership, and reuses existing
    graph/filter tasks. Category 5 remains open until real-data and final graph acceptance.
    The [regional-01 assessment](core-graph-regional-01-assessment.md) confirms clean
    execution and coordinate audits, but the one-line correction changes none of
    the cohort's surviving coverage or chromosome memberships. It is not accepted
    as the solution to missing arms. Distinguish the mismatched-field bug candidate
    from explicit interval-length policy; inspect rejection gates before further changes.
    The [observation-only gate batch](core-graph-regional-gates-test.md) now records
    pre/post-overlap evidence, individual candidate gates and the final coverage
    decision for both grouping variants. It requires byte-identical assignment
    outputs with/without observation and retains the existing graph/filter cache.
    [Gate-run acceptance](core-graph-gates-01-assessment.md) passed: the restored
    chr7 arms form 32.8/40.2 Mb runs at ~95% density but fail the explicit >1 Mb
    individual-interval gate (longest 49.8/102.6 kb). Corrected grouping works, but
    this independent policy keeps the final result unchanged. Compare supported
    chunk settings and all-cohort side effects before changing production defaults.
    The [supported chunk-setting batch](core-graph-regional-chunks-test.md) compares
    1 Mb / 100 kb / 25 kb under both grouping versions and all three filtering policies.
    Six 1 Mb baselines remain cached; 12 new assignment jobs and a report measure arm
    recovery, control losses, new small chromosome assignments and coordinate ownership.
    [Chunk-run assessment](core-graph-chunks-01-assessment.md) selects guarded
    filtering + corrected grouping + 25 kb as the production candidate: both missing
    arms recover, 98.74 Mb mapped gain, no source mapped losses or ordinary-control
    chromosome changes. The coupled softclip behavior adds 37.81 Mb unreported source
    sequence, so production integration must explicitly account for source ends and
    pass real FULL/CLIP/node audits before category acceptance or diagnostic cleanup.
    The selected policy is now integrated into normal production construction,
    independently of the removable diagnostics. Native tool regression, workflow,
    source-coordinate and ordinary resume checks pass locally. The
    [production rebuild and acceptance run](core-graph-production-01-test.md) now
    checks actual chromosome-arm representation, excluded source ends and distinct
    graph-node accounting. Category 5 remains open pending Crest acceptance.
    [Production-01 failed](core-graph-production-01-assessment.md) at a pinned
    Cactus fragment FASTA/PAF naming mismatch, before completed graph export.
    The production adapter now reconciles fragment identifiers without changing
    coordinates and probes the installed extraction/sanitizing tools before
    graph construction. Production-02 must pass the existing recovery/accounting
    acceptance checks; no assembly or chimera changes are indicated.
- [ ] **6. Full Markdown report** — integrate current-run artifacts, decisions, plots and statuses.
- [ ] **7. Operational cleanup** — dependency/resource/cache/output cleanup and final integration acceptance.
- [ ] **8. PSMC** — implement and accept both agreed input modes and reporting.
- [ ] **9. BlobToolKit evidence and cleaning** — evidence first, separately validated cleaning.
- [ ] **10. Graph self-mapping QC** — bounded supported mapping, metrics and report integration.
- [ ] **11. Private/accessory sequence evidence** — repair/replace and decide retained scope after requested use.
- [ ] **12. Graph SV interpretation and supplementary SyRI** — method validation and integrated reporting.
- [ ] **13. Additional graph exports/figures** — inventory and disposition, then only justified agreed additions.

Category 1 does not wait for optional post-assembly comparison plots or the consolidated
report. Their runtime/routing/dependency acceptance belongs to categories 2/7 and their
report integration to category 6. If later work exposes a concrete defect affecting an
accepted category, record and address it rather than silently retaining an invalid checkmark.

**2026-10-25 continuation:** completed pangenome replay/trace diagnostics are removed
from the active pipeline. Permanent production ownership/coordinate/native checks
remain. The [tool-change record](../docs/pangenome-tool-changes.md) distinguishes
implementation defects from policy choices and documents biological rationale,
validation and limitations. The [input/failure/cleanup checkpoint](core-input-failure-cleanup-checkpoint.md)
records the expanded bounded matrix and the specific remaining short-read,
taxonomy/reference, operational and reporting checks. Categories 2/3/6/7 remain open;
PSMC remains the first extension after their core acceptance gate.

The subsequent [organelle/taxonomy contract pass](organelle-failure-contracts-2026-10-25.md)
repairs swallowed GetOrganelle/MITOS failures, retry graph discovery, incomplete
versus circular naming, missing context joins and invalid taxonomy records. It adds
enabled-route fixtures and explicit recovery states; native short-read pilot and
multi-record annotation support remain unaccepted, so categories 2/3/7 stay open.

### 2026-10-16 input/reference/graph check batch

- Shared final linear-reference selection now runs independently of graph admission and
  graph enablement. Singleton and collapsed/fragmented assemblies can have a linear reference
  without becoming graph eligible. An eligible graph uses the same selected reference.
- Optional `--linear_reference_ids` takes comma-separated final assembly IDs, at most one per
  species. Unknown, missing or duplicate species overrides fail. An override unsuitable for
  the graph withholds that graph rather than silently choosing a different coordinate frame.
- `references/<taxid>/reference.json` records the final FASTA SHA256, sequence dictionary,
  source assembly/path, selection reason and mapping compatibility requirement. The real
  FASTA and FAI are checked together. This does not implement PSMC or certify biology.
- Final haplotype/ploidy metadata mismatches now fail eligibility. Graph input auditing no
  longer overwrites an existing member taxid with the cohort taxid, which could mask a
  mismatched species record. Both input and build reports reach the report channel.
- Explicit `mito_tool=none` now bypasses organelle assembly/reference lookup and emits empty
  bait lists so nuclear refinement continues for both HiFi and short-read routes.
- Use ordinary Nextflow caching and `-resume`. No extra persistent state, job-store recovery
  mechanism or protection against cache deletion was added. Successful cached tasks are
  reused; failed/missing tasks rerun through normal Nextflow behavior.
- The review archive now includes reference/eligibility records and compact graph reports,
  summary tables and figures (under 20 MB each), excluding graphs, mappings, VCFs and indexes.
- Acceptance remains scoped: synthetic production-wiring tests do not certify a real Cactus
  graph or every supported input route. Current runtime failures still stop the controller
  after configured retries. Independent-species continuation, short-read tool-failure versus
  no-recovery reporting, reference-lookup fallback and complete real short-read integration
  remain open under categories 2/3. Do not mark those categories complete from these checks.

The next real run enables construction with `--run_pangenome true
--pangenome_validate_only false`, retaining the reviewed cut file, current input sheets,
normal work/cache paths and all-stage QC. Its graph acceptance remains category 5.

Local acceptance: 47 focused Python tests passed (eligibility, reference identity, graph
input/output/variant audits, sharing/regional calculations and FASTQ subset preparation).
`python3 tests/run_core_contract_checks.py` passed ten Nextflow 23.10.1 scenarios, including
eight expected failures. External Cactus/MitoHiFi executables are tiny doubles; production
workflow modules and auditors execute normally. These tests cover routing/validation,
not the biological quality of real graph or organelle outputs. Replaying the uploaded final
name-map lengths preserved all ten graph identities, reference chromosomes and the existing
Sde-CMat_203_hap1 reference. Logs are under `comparisons/core-check-*.log`.

After syncing, from the Crest project directory, with no other controller using its resume session:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-graph-01 \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json" \
  --chimera_break "$PWD/data/chimera_reviews/Sde-CTlk_104.v7.reviewed.tsv" \
  --run_pangenome true --pangenome_validate_only false
```

This uses normal pipeline resources, work and resume paths, all-stage QC and symlink
publication. It retains the manual cuts. No new assembly parameter settings are introduced.
Return `comparisons/core-graph-01.review.tar.gz` for real graph/output acceptance.

### 2026-10-19 graph representation follow-up

The supplied CTlk name maps match the accepted final maps. Both haplotypes contain
chr4_1+chr14_1 and chr7_1+chr12_1, plus separately named chr4/chr7 pieces. Regional
statistics group whole scaffold names; they do not paint the components of a composite
onto separate chromosomes. Thus pure chr4/chr7 rows cover only separately named pieces.

An additional discrepancy was found: CTlk hap1's chr4_1+chr14_1 is absent under that
sequence name from the accepted CLIP path-group ledger, although it exists in the final
assembly. Hap2's counterpart and both chr7/chr12 composites are reported. This requires
FULL/CLIP path and clipping evidence; it must not be described as biological absence.

Implemented for the next regular resume:

- A small regional-context process compares final chromosome-scale selections with the
  existing regional measurements. Main figures include chromosome-scale composites and
  distinguish composite-only (`C`), separate pieces plus composite (`+`) and an expected
  composite missing under its name from CLIP (`!`, or `+!` with separate pieces).
- Composite composition is shown per haplotype/scaffold, including unavailable expected
  entries without fabricating percentages. Other composite labels remain in supplementary
  tables; no graph sequence or global sharing measurement is removed or reassigned.
- A separate post-build QC process reads the existing Cactus stats bundle, summarizes
  retained path lengths, unions reported clip intervals within tool path names, compares
  FULL/CLIP graph sizes and reports per-sequence FULL W-coordinate spans versus CLIP path
  lengths. Unmapped tool names are explicit; residual source lengths are not called clipping.
- Correct Cactus package/commit, Toil and verification-container provenance is captured
  separately. The original construction version record is retained; its Toil value is no
  longer propagated as the Cactus release in the workflow version channel.
- Statistics requests 128 GB; node coverage/similarity request 96 GB, matching observed
  memory with headroom. Review archives include `.nwk`, compact `.stats.tgz`, representation
  summaries and the regional context/presence tables.
- Cactus construction scripts, input assemblies, retained cuts and assembly parameters
  are unchanged. Ordinary Nextflow resume was verified with the real construction/QC
  modules and external tool doubles; the construction tasks were cached.

Validation: 60 focused Python tests passed. Ten Nextflow scenarios passed with the new
post-build/context processes; a further normal-resume check passed. Updated figures were
rendered and visually inspected from the actual run tables in `comparisons/core-graph-02-preview`.
The real Cactus clipping bundle is on Crest and still needs the new production QC step.
Category 5 remains open pending that evidence, especially the missing CTlk hap1 composite.

After syncing, from the Crest project directory:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-graph-02 \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json" \
  --chimera_break "$PWD/data/chimera_reviews/Sde-CTlk_104.v7.reviewed.tsv" \
  --run_pangenome true --pangenome_validate_only false
```

Return `comparisons/core-graph-02.review.tar.gz`. New results include
`results/pangenome/373251/representation/clipping_report.md` and
`results/pangenome/373251/regional/figures/regional_figures.md`.

### Existing post-assembly plots: acceptance still needed

Dotplots/riparian, cohort QUAST, final telomere visualization and coverage books have existing
implementations. Expect reuse, with bounded checks of final-input identities, pairing/species
scope, empty/disabled routes, dependencies and plot content before full report integration.
They were disabled in the recent accepted assembly/QC runs. Current all-pairs selection
builds pairs across its input stream; verify species-aware behavior for general multi-species
use. SETUP_PAFR still performs an unpinned runtime GitHub install and needs dependency cleanup.
For the current ten assemblies, all-pairs mode schedules 45 alignments; within-sample mode
schedules five. Their scope/cost should be deliberate when enabling the branch.

## Current accepted checkpoint

- General cohort chimera evidence, explicit manual cuts and downstream finishing have been
  exercised on the real ten-haplotype cohort. Automatic cutting remains explicitly deferred.
- Chromosome-only telomere extension, record-order restoration and all-stage QC passed the
  latest uploaded real run: 116 completed tasks, 661 cached, no failures/retries. Both reviewed
  CTlk hap1 cuts remain verified. Purging/biological optimization is deferred by the user.
- Final Hi-C mapping/filtering/metrics/contact maps have now been decoupled from
  run_post_assembly and tested locally in final_only/all_stages, QC-off and plots-off cases.
  Real Crest acceptance passed in run 1535116. Local wiring tests additionally cover
  final_only/QC-off/plots-off cases using tool stubs.
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

1. Completed 2026-10-16: final Hi-C QC plus validation-only pangenome preflight/input
   audit (run 1535116). Category 1 is checked off. No graph or PSMC was started.
   The review archive includes small graph input audits and final Hi-C plot/statistic products.
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


## Reporting follow-up — HiFi alternate placements (2026-10-22)

Evidence figures must show apparent spanning molecules versus credible spans after
checking same-haplotype alternatives, with ambiguous/rejected reads visually distinct.
Show unique-arm bridge counts and read-length feasibility at each assessed site.
Provide a compact read-placement view with alternative scaffold/locus, alignment
quality and clipping; distinguish same-haplotype competition from other-haplotype
similarity. Preserve per-library Hi-C context and applied/retained decision labels.
This is deferred to the reporting pass; the underlying placement tables already exist.

## Operational resume accepted (2026-10-25)

`core-operational-01` passed: two completed tasks (native taxonomy guard and run
summary), 855 cached, no failures/retries. All 111 packaged pangenome files and
160 common reference/QC/name-map files match accepted production 02 byte-for-byte.
The graph and all five HiFi organelle assemblies stayed cached. Compact organelle
outputs were missing from the review archive; packaging now includes their status,
annotations, logs and figures. Categories 2/3/7 remain open for the native short-read
pilot and remaining contracts; this all-HiFi run does not validate that branch.
See [run assessment](core-operational-01-assessment.md).

## Native short-read pilot correction (2026-10-26)

CMat_061's reads-based organelle assembly completed with six fragments; the optional
graph-pruning retry exposed a legitimate empty-target early return without the
usual completion marker. The wrapper now verifies that output state explicitly.
The MITOS adapter now preserves/annotates every fragment instead of rejecting
multi-record input. Native acceptance remains pending the resumed pilot; this
supersedes the earlier single-record-only limitation, not the requirement for
native evidence. See [retry/fragment correction](organelle-native-retry-fix-2026-10-26.md).

### 2026-10-27 native pilot result

[Native pilot review](core-shortread-native-01-assessment.md): SPAdes, Redundans,
decontamination and finalization completed, with native fragmented organelle
annotation retained. The overall launcher still failed because empty harmonization
references suppressed eligibility/summary work. That empty-list collection and
a disabled-chimera staging collision are repaired and covered by two actual-wiring
integration cases. Resume acceptance remains pending; the 6.36 Mb nuclear test
output and 4,945 bp partial mitogenome do not establish biological completeness.

**Completion acceptance:** the replacement packet, run 1553560, exits 0 with three
new tasks and 15 cached, no failures/retries. Eligibility, shared reference and
summary now execute correctly for the unharmonized shotgun-only input. The native
short-read execution gate is accepted; no further identical pilot rerun is needed.
See [completion assessment](core-shortread-native-01-completion-assessment.md).
Categories 2/3/7 retain only their other documented outstanding contracts; this
acceptance does not certify the biological completeness of the reduced test output.

### 2026-10-28 identity/reference chunk

Scientific-name-only input handling is now wired into main.nf; unknown/ambiguous
names and missing identities fail explicitly. Reference retrieval validates its
outputs and distinguishes deterministic invalid outputs/dependencies from bounded
external retries. Seventeen local workflow cases, a full-main input-entry check
and the native-check workflow with controlled tools passed. The small
[native identity/reference check](identity-reference-contracts-2026-10-28.md) is
ready for Crest; no assembly or graph rebuild is involved. Categories 2/3 remain
open pending that evidence and final checklist reconciliation; the broader retry
and dependency audit remains category 7.

**Native acceptance:** `core-identity-reference-01` exits 0 with four completed
tasks and no failures/retries. Scientific-name and numeric routes agree on 373251;
the reference is NC_009588.1 (16,617 bp). Categories 2 and 3 are now checked off
within their supported execution/contract scope. See
[native assessment and remaining gates](core-identity-reference-01-assessment.md).

## 2026-10-28 optional finalized comparisons

The production comparison branch has been repaired and locally validated; see
[contracts, fixes and native test](final-comparison-contracts-2026-10-28.md).
This is category 7 work, not reopening accepted input/organelle/graph categories.
The next bounded Crest check compares accepted final assemblies without rebuilding
them. Native comparison acceptance and the broader operational cleanup remain open;
unified Markdown reporting follows under category 6.

### Native comparison execution result

`core-final-comparisons-01` completed 13/13 tasks with no failures or retries.
Final input fingerprints and five within-sample pairs reconcile. Sparse retained
alignment evidence and presentation limitations prevent complete visualization
acceptance; see [review and bounded follow-up](core-final-comparisons-01-assessment.md).
No assembly or graph rebuild is indicated.

### R synteny restoration (2026-10-28)

The optional pairwise figures again use pafr/gggenomes with ggplot2 and the
established minimap2 mapping-chain mode. The Matplotlib renderer is retired.
Raw/filtered mapping retention and interval-union coverage are now measured;
R rendering, workflow cases and cache reuse passed locally. Native recovery of
expected mapping coverage remains open. See [restoration and one-pair Crest command](final-comparisons-r-restoration-2026-10-28.md).
The user prefers R/ggplot2 for plotting throughout the pipeline. Remaining active
Python figure renderers are migration items for the category 6 reporting pass.
