# gcl_genome_assembly — Multi-Input Refactor Plan

**Purpose.** Control document for expanding `gcl_genome_assembly` from a HiFi+Hi-C diploid-only
pipeline to one that accepts heterogeneous sequencing inputs (HiFi, Hi-C, TellSeq linked reads,
short-read shotgun) in any valid combination, produces haploid or diploid assemblies, and runs
only the steps appropriate to the data + assembly type.

**How to use this doc.** Update the Status Dashboard and the per-phase checkboxes as work lands.
Record completed work in the Change Log with a date. When a decision changes, edit the relevant
architecture section and add a note under §7. Keep this committed in the repo so any session can
pick up where the last left off.

**Status legend:** `[TODO]` not started · `[WIP]` in progress · `[DONE]` complete · `[DEFERRED]` parked
Task checkboxes: `- [ ]` open · `- [x]` done.

---

## Status Dashboard

| Field | Value |
|---|---|
| Current phase | Phase 0 — design (this doc) |
| Phases complete | none |
| Last updated | 2026-06-29 |
| Production branch | `main` (S. delicatulus runs — DO NOT break) |
| Refactor branch | TBD (e.g. `feature/multi-input`) |

---

## 1. Goal & Scope

**Supported input modes (target):**

1. HiFi + Hi-C (current behavior — must remain identical)
2. HiFi only (no scaffolding; final = conditioned contigs)
3. HiFi + TellSeq (linked-read scaffolding; final = linked-read scaffolds unless Hi-C also present)
4. Short-read only (SPAdes contig assembly; scaffolded only if Hi-C and/or TellSeq present)
5. Any of the above ± Hi-C ± TellSeq (combinations resolved per-sample)

**Cross-cutting:**

- Haploid **or** diploid output, driven by `meta`, not hardcoded.
- Per-sample heterogeneity within a single run (sample A HiFi+Hi-C diploid, sample B short-read haploid).
- 2+ rounds of YaHS Hi-C scaffolding retained and made easier to extend (round count is **global**).
- Downstream tail (`FINALIZE_ASSEMBLY` → TIDK, snail, QUAST_FINAL, pairwise, report, coverage book)
  reused unchanged once branches converge on the canonical assembly channel.

**Explicitly out of scope (for now):** ONT long reads, trio/parental binning, pangenome graph
assembly. Leave meta fields extensible so these can be added later.

---

## 2. Supported Input Matrix

Cell = tool/behavior for that mode. `—` = step skipped.

| Input mode | Contig assembler | Dedup* | Mito | Correct/polish† | Scaffolding | Gap-fill‡ | Teloclip | Contact maps |
|---|---|---|---|---|---|---|---|---|
| HiFi + Hi-C | hifiasm | purge_dups | MitoHiFi | Inspector | YaHS ×N | TGSGapCloser | yes | yes |
| HiFi only | hifiasm | purge_dups | MitoHiFi | Inspector | — | — | yes | — |
| HiFi + TellSeq | hifiasm | purge_dups | MitoHiFi | Inspector | ARCS/LINKS | TGSGapCloser | yes | — |
| HiFi + TellSeq + Hi-C | hifiasm | purge_dups | MitoHiFi | Inspector | ARCS/LINKS → YaHS ×N | TGSGapCloser | yes | yes |
| Short-read only | SPAdes | redundans | MitoFinder | Pilon | — | — | — | — |
| Short-read + Hi-C | SPAdes | redundans | MitoFinder | Pilon | YaHS ×N | Sealer | — | yes |
| Short-read + TellSeq | SPAdes | redundans | MitoFinder | Pilon | ARCS/LINKS | Sealer | — | — |
| Short-read + TellSeq + Hi-C | SPAdes | redundans | MitoFinder | Pilon | ARCS/LINKS → YaHS ×N | Sealer | — | yes |

\* **Dedup is user-selectable per sample** (`purge_dups | redundans | none`) **regardless of input**;
cells show the per-assembler default. redundans is always run in **reduction mode only**.
† **Correct/polish is optional** (param-gated, like decontamination). Cells show the tool used when enabled.
‡ Gap-fill runs only if scaffolding produced gaps. Teloclip is long-read-only by mechanism → never runs without long reads.

---

## 3. Target Architecture & Key Decisions

### 3.1 Meta-map carrier (foundational)

Replace positional tuples (`sample_id, hifi_bam, hic_r1, hic_r2`) and (`sample_id, hap1, hap2`)
with `tuple(meta, files)`. Compute meta **once** in the parser; carry it everywhere; stop deriving
identity via `replaceAll(/_hap[12]$/, '')`.

```groovy
meta = [
  id:         'CLim_110_hap1',   // unique per assembly unit
  sample:     'CLim_110',
  haplotype:  'hap1',            // 'hap1' | 'hap2' | 'primary'
  n_hap:      2,                 // 1 = haploid/collapsed, 2 = diploid  (drives groupKey)

  // read evidence present (booleans, set by parser from which columns/cells are populated)
  hifi:       true,
  hic:        true,
  tellseq:    false,
  shortread:  false,
  long_reads: true,              // derived (hifi || future ONT) — gates teloclip / long-read gap-fill

  // assembly strategy (defaults derived from evidence; overridable per row)
  assembler:  'hifiasm',         // 'hifiasm' | 'spades' | (future: more)
  dedup:      'purge_dups',      // 'purge_dups' | 'redundans' | 'none'  (selectable, input-independent)
  mito_tool:  'mitohifi',        // 'mitohifi' | 'mitofinder' | 'none'

  // scaffolding plan (derived from evidence + params)
  hic_rounds: 2,                 // N rounds of YaHS (0 if no Hi-C); GLOBAL param value
  scaffolders:['linked','hic']   // ordered; subset of available given evidence
]
```

### 3.2 Ploidy / groupTuple fix

Every `.groupTuple(by: 0, size: 2)` becomes `.groupTuple()` keyed on `groupKey(meta.sample, meta.n_hap)`.
This lets haploid (group of 1) and diploid (group of 2) coexist. A literal `size: 2` will hang
forever on a haploid sample.

### 3.3 Routing rule

- **Global capability toggle → `if (params.x)`** (e.g. QC mode, round count).
- **Per-sample, data-driven divergence → `.branch{}` + `.mix()`** (e.g. some samples have Hi-C, some don't).

Both are used. The `.branch{}` arms route to different subworkflows, each emitting the same contract,
then `.mix()` reconverges before the next shared stage.

### 3.4 Output-contract discipline

Every assembler branch and every scaffolder branch emits **`tuple(meta, fasta)`** per assembly unit.
Nothing downstream inspects which branch ran — it reads `meta` flags. This is what lets the
post-`FINALIZE_ASSEMBLY` tail stay untouched.

### 3.5 Hi-C N-round scaffolding

Encapsulate one round — remap (`MAP_HIC_TO_*`) + refilter (`FILTER_HIC_BAM_*`) + `SCAFFOLD_HIC`
(YaHS) + that round's pairs/BAM checkpoint QC — into a single subworkflow `SCAFFOLD_HIC_ROUND`
taking `(meta, assembly, hic_reads, round_params)` and emitting `(meta, scaffolds)` plus QC.

`main.nf` then chains explicit calls (compile-time, not a runtime loop):

```
r1 = SCAFFOLD_HIC_ROUND(contigs, reads, params.yahs_round1)
r2 = SCAFFOLD_HIC_ROUND(r1.scaffolds, reads, params.yahs_round2)
// add round 3 = one more call
```

- Round count is a **global param** (`params.hic_scaffold_rounds`, default 2). **Decided.**
- Per-sample N (if ever needed later): gate each round's input `.filter { meta.hic_rounds >= k }`, mix early-stoppers back at the convergence point.
- Runtime-variable N via Nextflow preview recursion: **declined** on 23.10.1 (fragile, version-sensitive).
- Bonus: each round self-contains its checkpoint QC, so N rounds → N metric sets automatically.

### 3.6 Contig conditioning subworkflow (data-type dispatch)

New subworkflow `CONDITION_CONTIGS(meta, contigs, reads…)` run immediately after assembly,
**before scaffolding**. Branches internally on `meta.assembler` / read type and emits clean contigs
+ QC. This is where the long-read vs short-read divergence lives.

Order within each branch (conceptual): `dedup → (optional) correct/polish → mito-removal → (optional) decontam`

- **Each step is independently optional / param-gated**, the same way contig decontamination already is.
- **dedup** is a selectable method (`purge_dups | redundans | none`) chosen via `meta.dedup`,
  **independent of read type**. Default derived from assembler (hifiasm→purge_dups, spades→redundans),
  overridable per row. redundans = **reduction mode only**. purge_dups read source follows available
  reads (long preferred); efficacy on short-read assemblies to be **verified at Phase 4**.
- **correct/polish**: long-read → Inspector (`--skip_baseerror`; already gated by `inspector_run_on_contigs`).
  Short-read → Pilon (new optional gate). Structural misassembly breaking from short reads alone is weak;
  let structural fixes come from Hi-C / linked-read evidence in scaffolding.
- **mito removal**: long-read → MitoHiFi; short-read → MitoFinder. Both feed the existing
  alignment-and-filter `FILTER_MITO_CONTIGS` pattern, so that step stays identical. Gated on `meta.mito_tool != 'none'`.

### 3.7 Post-scaffolding finishing gates

Run near finalize, gated on evidence/results, not on assembler:

- **gap-fill**: only if scaffolding ran (gaps exist). Long-read scaffolds → TGSGapCloser; short-read scaffolds → ABySS-Sealer.
- **teloclip**: only if `meta.long_reads`. Independent of scaffolding.
- **contact maps / compartments / TADs**: only if `meta.hic` AND `params.make_final_contact_maps`.

### 3.8 QC read-source generalization

`BUILD_MERYL_DB` and `MAPPING_QC` currently hardwire `BAM_TO_FASTQ.out` (HiFi). Generalize to a
"qc_reads / kmer_reads" channel populated by the input layer: HiFi for long-read assemblies, PE
shotgun for short-read assemblies, with the matching mapper (minimap2 vs bwa-mem2). Merqury QV uses
whichever read set is the truth set.

### 3.9 QC mode toggle (cost control)

`params.qc_mode = 'all_stages' | 'final_only'`.

- **`all_stages`** (default): full QUAST/BUSCO/MERQURY/MAPPING battery at every intermediate stage
  (= current behavior; no change).
- **`final_only`**: QC battery runs only on the finalized assembly (post-`FINALIZE_ASSEMBLY`,
  the canonical converged channel); the ~10 intermediate `ASSEMBLY_QC_*` calls are skipped.
- The final assembly is **always** QC'd in both modes.

Implemented as `if`-gating around the intermediate `ASSEMBLY_QC_*` calls. Default preserves behavior,
so this lands in Phase 1.

---

## 4. Tool Slots by Assembly Type

`(revisit)` = chosen for now, re-evaluate at the noted phase.

| Role | Long-read (HiFi) | Short-read | Notes |
|---|---|---|---|
| Contig assembly | hifiasm | SPAdes | selector built pluggable; more assemblers (long & short) possible later |
| Haplotig/dup reduction | `purge_dups \| redundans \| none` | `purge_dups \| redundans \| none` | **selectable, input-independent**; redundans = reduction-only |
| Mito assembly + removal | MitoHiFi → FILTER_MITO_CONTIGS | MitoFinder → FILTER_MITO_CONTIGS | optional (`mito_tool`); removal pattern shared |
| Base correction / polish (optional) | Inspector (`--skip_baseerror`) | Pilon `(revisit P4)` | param-gated |
| Structural misassembly | Inspector | from Hi-C / linked-read evidence | no strong SR-only tool |
| Gap-fill (post-scaffold) | TGSGapCloser | ABySS-Sealer `(revisit P5)` | only if scaffolding ran |
| Teloclip | yes | n/a (needs long reads) | — |
| Merqury / mapping QC reads | HiFi | PE shotgun | follows assembler |

**Future / revisit:** redundans as a **scaffolder / gap-closer** option, used alongside linked-read +
Hi-C scaffolding — decide how to factor into the scaffolding process at **Phase 5**.

---

## 5. Phased Implementation Plan

Each phase: branch from the previous, validate on **one test sample**, then merge. Changing tuple
shapes invalidates all cache hashes — expect a full re-run per phase (don't rely on `-resume`
carrying across a phase boundary).

> **Phase 4 and 5 swapped** (vs original plan): short-read contig path comes first because shotgun
> short-read test data is available now; linked-read scaffolding comes last because no linked-read
> test data exists yet.

### Phase 1 — Meta-map refactor (behavior-identical) `[TODO]`

Goal: introduce the meta carrier, `groupKey`, and the QC toggle with **zero behavior change** on HiFi+Hi-C diploid.

- [ ] Rewrite `parse_sample_sheet.nf` as a validating parser emitting `tuple(meta, files)`; keep `Channel.fromList()` for per-sample cache independence; keep relative-path resolution.
- [ ] Define the canonical `meta` shape (§3.1) in one place (a helper or the parser).
- [ ] Thread `meta` through `BAM_TO_FASTQ`, `TRIM_HIC`, `HIFIASM`, mito filtering, purge, Inspector, scaffolding, gap-fill, teloclip, finalize — updating each module's input/output signature.
- [ ] Replace every `groupTuple(by: 0, size: 2)` with `groupKey(meta.sample, meta.n_hap)` grouping.
- [ ] Replace `replaceAll(/_hap[12]$/, '')` identity surgery with `meta` field reads.
- [ ] Add `params.qc_mode` (`all_stages` default | `final_only`); `if`-gate the intermediate `ASSEMBLY_QC_*` calls so default reproduces current behavior (§3.9).
- [ ] Validation: run one S. delicatulus sample; confirm outputs match the current `main` results.

### Phase 2 — Ploidy generalization `[TODO]`

Goal: make diploid a meta-driven `n_hap`; support haploid output.

- [ ] Drive hifiasm `--primary` and `n_hap=1` from `meta.ploidy`.
- [ ] Confirm all grouping/QC paths handle `n_hap=1` (no hang, no missing-pair errors).
- [ ] Gate hap1-vs-hap2 pairwise/riparian on `n_hap == 2`.
- [ ] Validation: run an existing sample forced haploid; one assembly out, full QC, no stalls.

### Phase 3 — Input flexibility + assembler selector + QC read source `[TODO]`

Goal: flexible sample sheet, contig-assembler selection scaffold, generalized QC reads. Enables HiFi-only.

- [ ] Wide sample sheet, **header-driven**: a missing column **or** an empty cell = no data of that type for that sample. Parser derives `hifi/hic/tellseq/shortread` flags and validates combinations (reject "Hi-C scaffolding requested but no Hi-C reads", "no read types", etc.). **(decided)**
- [ ] `CONTIG_ASSEMBLY` selector subworkflow emitting `(meta, fasta)`; hifiasm wired, short-read branch stubbed (filled in Phase 4); built to accept future assemblers.
- [ ] Generalize `BUILD_MERYL_DB` + `MAPPING_QC` to a qc_reads channel (§3.8).
- [ ] Gate `MITOHIFI` on `meta.hifi`.
- [ ] Validation: HiFi-only run (drop Hi-C columns) → conditioned contigs as final, no scaffolding, no contact maps, teloclip still runs.

### Phase 4 — Short-read contig path `[TODO]`  *(was Phase 5)*

Goal: SPAdes assembly + a real branched `CONDITION_CONTIGS`. Validated **short-read-only**.

- [ ] Wire SPAdes into the `CONTIG_ASSEMBLY` selector (Jason has scripts).
- [ ] Build `CONDITION_CONTIGS` subworkflow (§3.6) with internal branching:
  - long-read branch wraps existing purge_dups / Inspector / MitoHiFi-removal (each optional);
  - short-read branch = dedup (`redundans|purge_dups|none`) → optional Pilon → MitoFinder-removal.
- [ ] Make `dedup` method selectable via `meta.dedup`, input-independent; redundans reduction-only.
- [ ] **Verify** purge_dups behaves sensibly on a short-read assembly (read source + duplicate calls).
- [ ] Confirm long-read-only steps (Inspector, teloclip) are skipped for SR-only; merqury/mapping use PE reads.
- [ ] **Regression**: re-run the HiFi+Hi-C production sample through the new `CONDITION_CONTIGS`; outputs must still match `main`.
- [ ] Validation: short-read-only run → SPAdes contigs, conditioning, full QC, no scaffolding, no teloclip, no contact maps.
- [ ] **Note:** short-read+Hi-C and short-read+TellSeq combos validate in **Phase 5** (need scaffolding chain).

### Phase 5 — Scaffolding chain + linked reads `[TODO]`  *(was Phase 4)*

Goal: per-round Hi-C subworkflow, short-read gap closer, linked-read scaffolder, data-driven routing.

- [ ] Refactor Hi-C scaffolding into `SCAFFOLD_HIC_ROUND` (§3.5); express round1/round2 as chained calls; add `params.hic_scaffold_rounds` (global).
- [ ] `.branch{}` routing by evidence: linked → hic → none; `.mix()` to converge on `ch_final_assembly`.
- [ ] Apply finishing gates (§3.7): gap-fill (TGSGapCloser long-read / ABySS-Sealer short-read), teloclip, contact maps.
- [ ] Add linked-read scaffolding module (ARCS/LINKS lineage, optionally Tigmint) emitting `(meta, fasta)`. **Build to contract; live validation DEFERRED — no linked-read data.**
- [ ] **Revisit:** redundans as a scaffolder/gap-closer option alongside linked-read + Hi-C; decide how to factor in (see §4 Future).
- [ ] Validation: HiFi+Hi-C 3-round run works; short-read+Hi-C scaffolds and gap-closes; (HiFi+TellSeq / short-read+TellSeq deferred — contract-only).

---

## 6. Change Log (implemented)

_(Add dated entries as work lands. Newest first.)_

- _(none yet)_

---

## 7. Decisions & Open Questions

### Resolved (2026-06-29)

- **Sample-sheet shape:** wide CSV, header-driven; missing column **or** empty cell = no data of that type for that sample.
- **Scaffolding rounds:** global param (`params.hic_scaffold_rounds`, default 2).
- **Short-read assembler:** SPAdes (Jason has scripts). Selector built pluggable for future options.
- **Dedup:** selectable method (`purge_dups | redundans | none`), **independent of input type**; redundans reduction-only; default derived from assembler.
- **QC mode:** `params.qc_mode = all_stages | final_only`, default `all_stages` (current behavior).
- **Correct/polish:** optional (param-gated) for both long- and short-read paths.
- **Phase order:** Phases 4 and 5 swapped (short-read before linked-read scaffolding).
- **Conditioning structure:** two dispatch points — `CONDITION_CONTIGS` (pre-scaffold) + finishing gates (post-scaffold).

### Still open

- [ ] **purge_dups on short-read assemblies:** confirm read source + that duplicate calls are trustworthy (Phase 4).
- [ ] **Short-read polisher:** Pilon for now — re-evaluate at Phase 4.
- [ ] **Short-read gap closer:** ABySS-Sealer for now — re-evaluate at Phase 5.
- [ ] **redundans as scaffolder/gap-closer:** how to factor into the scaffolding chain alongside linked-read + Hi-C (Phase 5).
- [ ] **Additional assembler options** (more short-read; other long-read assemblers): pluggable selector, deferred.

---

## 8. Risks & Test Discipline

- **Cache invalidation:** any tuple-shape change rehashes every task; plan a full re-run on the test
  sample per phase. Branch in git; keep `main` runnable for production.
- **groupTuple hang:** the `size: 2 → groupKey` swap must be complete before any haploid run.
- **Conditioning refactor regression:** Phase 4 pulls long-read conditioning into a subworkflow —
  re-validate the HiFi+Hi-C production path after that change.
- **Heterogeneous-run correctness:** test at least one mixed-input sample sheet before declaring a phase done.
- **Decouple QC emission:** keep step subworkflows emitting `(assembly, qc_summary)` and let `main.nf`
  decide whether to thread the summary into `COMPILE_FINAL_QC` (preserves the current `.mix()` pattern).
- **One-tool-does-more surprises:** redundans (scaffolds/gap-closes — keep to reduction mode here),
  MitoFinder (annotates) — pin each slot to the behavior you actually want.
