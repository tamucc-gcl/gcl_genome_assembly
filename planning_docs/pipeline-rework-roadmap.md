# Pipeline rework: architecture and ordered implementation

Status: implementation roadmap; no source edits or execution performed. 23 September 2026.

This document supplements [the agreed biological and operational decisions](assembly-pangenome-agreed-plan.md) and [the static review](pipeline-static-review.md). Those decisions remain authoritative. The user runs all validation and biological jobs on a separate computer after syncing source through GitHub. The assistant edits and statically reviews source; it does not execute pipeline commands, tests, or analysis scripts.

## 1. Working sequence

Design the interfaces for the whole pipeline first. Establish workflow nesting, process aliases, parameter ownership, and output contracts before the first expensive validation baseline. Then implement upstream to downstream in independently reviewable batches. Each batch must remain runnable; do not deliver a partially wired replacement as a validation checkpoint.

The first structural migration may invalidate substantial existing cache. Avoid paying that cost repeatedly: complete foundational restructuring and necessary upstream fixes before asking for a full biological baseline. Existing cached work is useful only when Nextflow legitimately recognizes it; do not copy outputs into task directories or disguise changed inputs to force reuse.

An implementation checkpoint requires static review followed by the user's execution evidence before dependent scientific behavior is treated as validated. Independent downstream design may continue while the user runs a checkpoint.

## 2. Architecture

Keep `main.nf` as the entry point that connects major workflows and handles top-level run context. Move large inline closures, embedded report generation, and detailed branch orchestration into appropriate functions, scripts, or named workflows. Previous compilation failures are an architectural constraint, but their precise cause must be determined from the actual error if they recur; source line count alone is not a diagnostic.

Proposed responsibilities below are contracts, not a requirement to create one new file for every row. Reuse existing workflows where sensible and choose permanent names during batch 0.

| Boundary | Existing source anchors | Contract |
|---|---|---|
| Input preparation | `functions/parse_sample_sheet.nf`, `functions/meta.nf`, `nextflow.config` | Validated biological identities, read sets, effective configuration, and explicit skipped-row records |
| Read preparation and assembly | `workflows/contig_assembly.nf`, read-QC workflows, `modules/hifiasm.nf`, `modules/spades.nf` | Per-sample reads and metadata to assemblies with explicit representation and haplotype identity |
| Assembly refinement | Existing organelle/decontamination workflows; correction, purging, scaffolding and finishing modules | Versioned assembly stages, retained/removed sequence provenance, usable auxiliary data, and statuses |
| Assembly QC | `workflows/assembly_qc.nf`, `workflows/qc_phase.nf`, Hi-C QC workflows | Stage-specific QC artifacts under existing none/final_only/all_stages policy |
| Harmonization and chimera handling | `workflows/harmonize_scaffolds.nf`, `workflows/chimera.nf` | Species-keyed eligible cohorts, current coordinate maps, explicit bypasses, final assembly records |
| Pangenome construction | `workflows/pangenome.nf`, `modules/cactus_pangenome.nf`, manifest modules | Auditable cohort and reference selection; CLIP, GREF(CLIP), compact FULL support artifacts |
| Pangenome analyses | Existing pangenome modules, separated from construction orchestration | Named graph-role inputs to required tables, variants, QC, and optional analyses |
| Reporting | `workflows/reporting.nf`, report modules | Current-run artifact/status records to one GitHub Markdown report |

Keep individual tool execution in process modules. Workflows compose processes and apply documented routing; they should not grow into another monolithic entry point. Check all fully qualified configuration selectors when changing nesting or aliases.

### Shared contracts to settle before implementation

- Identity: species key, individual/sample ID, assembly ID, representation, haplotype ID, library ID and read-set ID are distinct. Never infer biological identity from a shortened filename.
- Reads: paired reads remain associated by explicit read-set identity and matching list order. The optional Hi-C table replaces the main-sheet Hi-C fields for listed samples. Independent libraries retain separate provenance.
- Assemblies: every record identifies its processing stage and coordinate version. Sequence-changing stages emit the information needed to invalidate or regenerate dependent mappings and coordinate products.
- Cohorts: expected eligible individuals are recorded before graph scheduling. Failed required processing cannot silently shrink a species cohort.
- Statuses: distinguish completed, skipped with reason, disabled, failed and unavailable. Missing files and synthetic placeholder successes are not interchangeable.
- Artifacts: named roles, owner, stage, path and availability feed reporting. Graph role is explicit rather than inferred from filename.
- Channels: declare input/output shape and cardinality for each workflow, including singleton and empty cases. Join by the complete biological key and assert unintended duplicate keys; avoid unkeyed cross-products and positional pairing.
- Task inputs: pass only stable metadata needed by that task. Keep growing report/provenance records separate from computational inputs so unrelated report changes do not invalidate assembly work.

## 3. Ordered implementation batches

### Batch 0 — Freeze architecture and migration contracts

Inspect `main.nf`, all active workflow includes, configuration selectors and parameter ownership. Create a process-to-workflow map and retained/replaced/removed inventory tied to the review findings. Record supported Nextflow/tool versions and the user's current execution profile. Establish permanent workflow names and narrow input/output contracts; move orchestration without deliberately changing scientific behavior in the same structural commit.

Prepare a lightweight status/artifact reporting contract now, even though final presentation is completed later. Do not postpone failure semantics and report inputs until the last batch.

Acceptance: all active includes, aliases, selectors and optional branches reconcile on static inspection; the user's compile/small-fixture check succeeds with the pinned runtime. Runtime evidence is required to verify compilation.

Cache expectation: broad invalidation is possible from changed workflow paths. Do not request an expensive full baseline solely to validate this structural step.

### Batch 1 — Inputs, identities and assembly routing

Replace naive delimited parsing; validate identity collisions, paths, paired inputs, overrides and skip reasons. Implement the optional multiple-Hi-C-library table. Resolve all read sets before route selection. Repair mixed HiFi/short-read routing and per-sample parameter propagation. Preserve Hi-C-only rows as announced assembly skips.

Acceptance: the user checks ordinary single-library inputs, the two-independent-library sample, HiFi-only, short-read-only, mixed-read, Hi-C-only, mixed-species and malformed/duplicate cases. Metadata reaches the correct route without cross-sample assignments.

Cache expectation: changed input contracts may affect all descendants. Complete this contract before validating expensive assembly tasks.

### Batch 2 — Read preparation and contig assembly

Implement per-library preparation and supported multi-file hifiasm input. Audit library-aware duplicate handling in its actual downstream consumer; do not add redundant processing. Repair effective assembler routing, representation labeling, thread/resource caps and pinned environments. Retain agreed assembly and purging defaults.

Acceptance: representative short-read and HiFi routes produce the intended assembly representation; primary/alternate is never mislabeled as hap1/hap2. Pair ordering, library provenance, filenames and sample cardinality are correct.

Cache expectation: input/read-preparation and assembly tasks may rerun. Once accepted, freeze their names, interfaces and scripts unless a correctness issue requires changing them.

### Batch 3 — Refinement, scaffolding and existing assembly QC

Repair organelle no-result behavior, conservative sequence separation, per-sample reduction controls and all conditional correction/scaffolding/finishing paths. Retain the agreed defaults and QC modes. Keep the nonfunctional BlobTools evidence branch disabled/unsupported; do not rebuild it here. Do not introduce annotation functionality.

Acceptance: skipped optional steps pass through the appropriate assembly; auxiliary recovery absence does not block nuclear output; failures remain visible. Multiple Hi-C libraries contribute correctly to one sample's downstream work. Every sequence-changing step has a consistent downstream coordinate/mapping contract. QC mode changes affect only their intended consumers.

Cache expectation: valid contig assembly should resume; modified refinement and its descendants rerun. Investigate unexpected upstream reruns before expanding the validation workload.

### Batch 4 — Harmonization, chimera handling and graph eligibility

Fix species-keyed matching, singleton/empty bypasses, current-coordinate cuts and preservation of nonparticipating routes. Evaluate chromosome-scale voting suitability from assembly evidence. Separate voting eligibility from pangenome membership. Enforce species reference requirements, minimum eligible individuals and required-failure cohort policy.

Acceptance: mixed species never share evidence accidentally; a single sample completes its assembly/report path; fragmented phased assemblies can be nonvoters without automatic exclusion; disabled chimera cutting preserves detection behavior; explicit cutting uses current coordinates.

Cache expectation: assembly/refinement should resume; affected harmonization, graph inputs and descendants rerun.

### Batch 5 — Pangenome construction and graph roles

Refactor construction separately from analysis. Establish CLIP as the biological graph, GREF(CLIP) as the default variant coordinate representation, and FULL as compact QC support. Repair species-specific index matching and biological path identities. Audit retry behavior and preservation of valid expensive construction work. Pin the chosen working toolchain rather than combining an untested upgrade with the refactor.

Acceptance: graph manifests reconcile to the intended cohort; synthetic gref paths do not increase biological denominators; no cross-species index/reference pairing occurs; construction artifacts satisfy downstream contracts. Previous-run comparisons use matched biological cohorts and corresponding events rather than raw VCF counts alone.

Cache expectation: graph construction can require a long rerun. Stabilize its required output contract before beginning that run. Later analysis-only changes should not modify this process merely for convenience.

### Batch 6 — Standard analyses and targeted method decisions

Verify Panacus exports against required per-haplotype, per-chromosome and chromosome-by-haplotype base-pair/percentage outputs. Replace overlapping bespoke summaries with published tools plus minimal necessary aggregation. Implement configurable sharing categories and one documented matrix for haplotype comparisons. Repair VCF tier routing, distance/missing-data treatment and retained temporary private-evidence calculations. Remove agreed unwanted branches once their required functions have replacements.

Resolve the INV/DUP/translocation method milestone here with user-run bounded experiments and published-tool comparisons. Keep this capability required without prematurely promising that one tool or gref alone supplies it. Self-mapping remains optional; individual population analysis remains out of scope.

Acceptance: tables reconcile to explicit denominators; chromosome accounting handles unplaced sequence; graph-sharing metrics exclude synthetic paths; variant coordinates/FASTA/sample identities agree; ordination labels match methods. Structural-variant replacement has evidence on known events and relevant assembly limitations before retiring the old implementation.

Cache expectation: accepted graph construction should resume. Rerun changed analyses and their reports, except when experiments reveal a necessary upstream correction.

### Batch 7 — Final reporting, publication and maintenance cleanup

Build one detailed GitHub Markdown report with collapsed analysis sections, concise visible run status and relative artifact links. Publish directly to the configured results directory. Use current-run manifests instead of directory globbing; introduce no automatic directory management. Remove retired code, obsolete parameters and dependencies only after the retained paths are validated.

Acceptance: single short-read, multiple same-species phased assemblies and mixed-species/mixed-input runs all produce meaningful reports. Disabled/skipped/failed states are clear, links resolve, and old outputs cannot masquerade as current results. Static inventory finds no active references to removed branches.

Cache expectation: presentation changes should primarily rerun reporting. A new computational requirement discovered here must be justified explicitly rather than quietly widening the rerun.

## 4. Validation and handoff for every batch

Each handoff records the source commit, changed files, intended behavior, static checks performed, user-run validation cases, expected cache reuse, expected reruns, acceptance criteria and known limitations. The user returns the command/configuration, Nextflow version, session identity, trace and relevant error/task logs as needed. No runtime success is claimed from static review.

Use small valid fixtures to exercise routing and known failure cases before representative biological inputs. Tiny fixtures cannot establish biological assembly quality. Check optional branches only where relevant to the batch, then perform final integration coverage across the agreed supported input combinations. Do not rerun every expensive case for an unrelated report edit.

For a failed checkpoint, repair the earliest responsible stage and rerun the bounded case. Do not advance dependent behavior as accepted or silently alter the biological defaults to obtain a pass. Preserve a known source commit for diagnosis; rollback does not guarantee cache reuse.

## 5. Cache and reproducibility rules

Nextflow's task identity depends on task/workflow names, inputs, scripts and execution environments. Resume also requires the task cache and surviving work outputs. The user should retain both on the execution computer and resume the intended session. Stable names do not guarantee reuse when scientific inputs change. See [Nextflow caching documentation](https://docs.seqera.io/nextflow/cache-and-resume).

Design choices for this rework:

- Freeze workflow nesting and aliases early; audit resource selectors at the same time.
- Keep input ordering deterministic, species/sample joins explicit and metadata task-specific.
- Avoid modifying staged inputs in place or embedding changing timestamps/report metadata in computational scripts.
- Select pinned tool/database versions before validating the stage that uses them.
- Keep report assembly independent of expensive computational modules; changing a report should not require editing their output declarations repeatedly.
- Preserve correctness over cache reuse. Never weaken cache validation to hide changed inputs.
- Maintain an expected-versus-observed cache table during user-run checkpoints. Investigate unexpected misses using the actual runtime's supported diagnostics.

DSL2 workflow boundaries follow [Nextflow's named workflow and input/output model](https://docs.seqera.io/nextflow/workflow). Use syntax supported by the pinned execution version rather than assuming every feature in current documentation is available.

## 6. Immediate next implementation step

Begin batch 0 with a source-level ownership map and precise channel contracts. Before changing pipeline source, read applicable repository instructions, inspect working-tree changes, and preserve unrelated user edits. Source write access may need a narrowly scoped filesystem permission because this repository is outside the current writable workspace. This roadmap does not request execution on the assistant's computer or a new biological run yet.
