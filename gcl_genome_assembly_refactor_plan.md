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
| Current phase | **Phase 2 (haploid) — scoping.** Fork + re-pairing already `n_hap`-aware; 3 design decisions (D1–D3, §7) pending before edits |
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
]
// NOTE (caching): scaffolding round count + scaffolder ordering are intentionally NOT in meta.
// meta is part of every task hash; keeping volatile knobs out means tuning
// params.hic_scaffold_rounds won't re-run upstream assembly/QC. Round count is read from
// params.hic_scaffold_rounds at point of use; scaffolder ordering is derived from
// meta.hic / meta.tellseq in the Phase 5 routing.
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

### Phase 1 — Meta-map refactor (behavior-identical) `[WIP]`

Goal: introduce the meta carrier, `groupKey`, and the QC toggle with **zero behavior change** on HiFi+Hi-C diploid.

- [x] Rewrite `parse_sample_sheet.nf` as a validating parser emitting `tuple(meta, files)`; keep `Channel.fromList()` for per-sample cache independence; keep relative-path resolution.
- [x] Define the canonical `meta` shape (§3.1) in one place (a helper or the parser).
- [x] **Caching:** exclude volatile fields (`hic_rounds`, `scaffolders`) from `meta` so tuning `params.hic_scaffold_rounds` does not re-hash/re-run upstream tasks under `-resume`.
- [ ] Thread `meta` through all module signatures, by pass:
  - [x] **Spine** — `bam_to_fastq`, `trim_hic`, `hifiasm` + the single haplotype fork (`ch_contigs`).
  - [x] **Mito** — `mitohifi`, `filter_mito_contigs` (collapsed HAP1/HAP2 → one process), `mito_circular_map` + STEP 4a rewire.
  - [x] **Purge + Inspector** — `purge_dups`, `correct_misassemblies` (contig) + collapse rewire onto `ch_mito_filtered`.
  - [x] **Decontam** — `decontaminate_assembly` subworkflow + 3 FCS modules + rewire (single per-hap `(meta, fasta)` channel).
  - [ ] **Hi-C map/scaffold spine** (2-round retained; `SCAFFOLD_HIC_ROUND` consolidation deferred to Phase 5):
    - [x] modules (`map_hic_to_assembly`, `filter_hic_bam`, `scaffold_hic`, `hic_mapping_metrics`) + contig-stage wiring (STEP 6–8) + contig checkpoints (`contig_raw_map`, `contig_filtered`, `scaffold_space`)
    - [x] scaffold-round2 (Inspector/decontam on scaffolds, remap/refilter, round 2) + round-2 checkpoints (`scaffold_round2_raw_map`, `scaffold_round2_filtered`, `scaffold_round2_space`)
  - [x] **Gap-fill / teloclip / finalize** — `gap_filling`, `teloclip`, `finalize_assembly` + rewire → `ch_finalized_assembly` (per-hap canonical channel feeding QC).
  - [x] **Final contact-map** — `contact_map` (+ `hic_compartments`/`hic_tads`) + STEP 14 rewire: `MAP_HIC_TO_FINAL`/`FILTER_HIC_BAM_FINAL` + `final_raw_map`/`final_filtered` checkpoints + `CONTACT_MAP_FINAL`. Consumes `ch_finalized_assembly`.
  - [ ] **QC pass** (final threading block):
    - [x] `assembly_qc` subworkflow + modules — `build_meryl_db`/`busco`/`mapping_qc` threaded (per-hap `meta`); `quast`/`merqury`/`combine_assembly_qc` unchanged (sample-string keyed); internal `groupKey` re-pairing
    - [x] 11 `ASSEMBLY_QC_*` call sites → pass per-hap channel (dropped the `groupTuple` prep blocks) + `params.qc_mode` gate on intermediates
    - [x] raw-read QC: `hic_qc`, `hifi_qc` threaded + call-sites (`HIC_QC_RAW` reads-map fix, `HIFI_QC` → `.out.fastq`); `hic_qc_from_bam`/`hic_scaffold_qc` are dead includes (see §9)
    - [x] downstream QC/report:
      - [x] STEP 14 viz tail — pairwise/riparian, `QUAST_FINAL`, TIDK call-sites (via `ch_final_by_id` string-id view; modules unchanged) ✅ applied
      - [x] post-QC compilation — metrics collection, `COMPILE_FINAL_QC`, `ASSEMBLY_REPORT`, `SNAIL_PLOT_FINAL`, `COVERAGE_BOOK`, manifest, `SUMMARY_REPORT` ✅
- [x] Replace every `groupTuple(by: 0, size: 2)` with `groupKey(meta.sample, meta.n_hap)` grouping. *(Done across all active code: the call-site re-pairs were deleted, `ASSEMBLY_QC` uses `groupKey`. No `size: 2` remains in `main.nf`.)*
- [ ] Replace `replaceAll(/_hap[12]$/, '')` identity surgery with `meta` field reads. *(Assembly-QC + spine done; 2 remain in the pairwise/viz tail — next sub-pass.)*
- [x] Add `params.qc_mode` (`all_stages` default | `final_only`); `if`-gate the intermediate `ASSEMBLY_QC_*` calls so default reproduces current behavior (§3.9).
- [ ] Validation: run one S. delicatulus sample; confirm outputs match the current `main` results.

### Phase 2 — Ploidy generalization `[IN PROGRESS — scoping]`

Goal: make diploid a meta-driven `n_hap`; support haploid (`n_hap==1`, `primary`) output.

**Already haploid-ready (verified in merged code):** `forkHaplotypeMeta` returns one `primary` meta for `n_hap==1`; `buildMeta` sets `n_hap=1` for haploid ploidy; the `assembly_qc.nf` `groupKey` re-pairing groups a variable hap count; pairwise `within_sample` `replaceAll(/_hap[12]$/)` won't match `primary` → no spurious hap-pair.

**Two-haplotype assumption to remove (dependency order):**
- [ ] **Fork call site** (`main.nf` ~566–571): `flatMap { meta, hap1_fa, hap2_fa -> [tuple(hmetas[0],hap1_fa), tuple(hmetas[1],hap2_fa)] }` — for `n_hap==1`, `hmetas[1]` is null. Zip `forkHaplotypeMeta(meta)` against the fasta(s) generically; haploid emits one `(primary, <primary fasta>)`. (Depends on D1: which fasta.)
- [ ] **`quast.nf` (QUAST)**: `input: tuple(sample_id, hap1_fasta, hap2_fasta)` + fixed `.hap1,.hap2` labels → variable `path(fastas)` + `val(labels)` (copy the `QUAST_FINAL` pattern already in this file).
- [ ] **`merqury.nf` (MERQURY)**: `merqury.sh db hap1 hap2 prefix` → build asm-arg list dynamically (`merqury.sh` takes 1 or 2 assemblies).
- [ ] **`assembly_qc.nf`**: replace `fastas[0], fastas[1]` (QUAST + MERQURY prep, lines 65/72) with the ordered fasta LIST + a labels list; feed the refactored modules.
- [ ] **`combine_individual_assembly_qc.R`**: `str_extract(..., 'hap[0-9]')` + `pivot_wider(names_from=haplotype)` → also recognize `primary`, tolerate a single hap. (Depends on D2/D3.)
- [ ] **`compile_qc.R`**: `pivot_longer(cols=c(hap1,hap2))` (×5) + `c(hap1,hap2,both)` → tolerate samples with only one hap. (Depends on D3.)
- [ ] Gate hap1-vs-hap2 pairwise/riparian explicitly on `n_hap == 2` (currently implicit via `replaceAll`).
- [ ] Validation: run a sample forced `ploidy=1`; one assembly out, full QC, no stalls; diploid parity unaffected.

**Design decisions (D1–D3) — see §7.**

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

- **2026-07-05 — Phase 2 stub validation PASSED (both paths, one run).** Local `-stub` (`stub_local.config`, `--decon_make_blobtools_evidence false`) on a 2-sample sheet — `Sde-CMat_203_dip` (ploidy 2) + `Sde-CMat_203_hap` (ploidy 1), same reads via distinct Hi-C symlinks — ran **192 processes to SUCCESS**. Confirms end-to-end: hifiasm fork emits 2 (diploid) / 1 (haploid) → per-hap steps `3 of 3` (dip_hap1+dip_hap2+hap_primary); QUAST/MERQURY 2-arity (dip) + 1-arity (hap); `groupTuple`/`COMBINE_ASSEMBLY_QC` `2 of 2` per stage; **pairwise/riparian ran only `dip_hap1_vs_dip_hap2`, no haploid pair (Pass 6 confirmed)**; COMPILE_FINAL_QC/SUMMARY_REPORT completed. Caveat: stub `touch`es outputs, so the R haploid logic (combine/compile/report) is NOT yet exercised — that's the real run's job. Two config-parse/wiring issues cleared en route (`params.outdir` default; `--mitohifi_species` required at build time). Next: real validation run.

- **2026-07-03 — Phase 2, pass 6 + 5b (haploid: within-sample gate + summary report).** **Pass 6** (`main.nf` ~1046): `within_sample` pairwise grouping now strips `_(hap[12]|primary)$` and documents that a haploid sample (one `primary` assembly) is intentionally excluded by the existing `haps.size()==2` check — so no within-sample dotplot/riparian is generated for haploid (already true implicitly via size==2; now explicit + robust). **Pass 5b** (`generate_summary_report.R`, 5 edits): sample-id derivation (L90), Section-1 assembly-files table (L163-170), and the within/cross dotplot count (L433-436) now strip/​map `primary`→`hap1`; Section-1 and Section-2 per-hap tables NA-fill a missing `hap2` column so a haploid row renders `value / —`; Section-2 diploid metrics (`qv`, `kmer_completeness`) fall back to the `hap1` value when there is no `both` line (haploid MERQURY single-asm); Section-3 visual loop resolves each sample's haplotype ids from the manifest (diploid → hap1/hap2, haploid → `primary` in slot 1, slot 2 = —) instead of hardcoding `paste0(sid,"_hap1"/"_hap2")`. Sections 4–7, `final_detail`, and every `across(c(hap1,hap2))` / `pivot_longer(c(hap1,hap2,both))` need no change — `compile_qc.R` guarantees those columns and the code uses tolerant patterns (`matches()`, `values_drop_na`). Diploid report output unchanged. **Phase 2 code-complete pending validation.**

- **2026-07-03 — Phase 2, pass 5 (QC metric R layer).** `combine_individual_assembly_qc.R`: added `norm_hap()` (hap1/hap2 → hap1/hap2; `primary` → `hap1` per D3; both/Both → both), applied at all 5 haplotype-extraction points (QUAST/BUSCO/MAPPING/merqury-qv/completeness), and now guarantees a uniform `hap1/hap2/both` output schema (NA-fills missing cols + fixed column order) so `compile_qc.R` is ploidy-agnostic. `compile_qc.R`: BAM + pairs `haplotype_id` parsing (lines ~106/113) now strip `_(hap[12]|primary)$` and map `primary`→`hap1`; plots/CSV unchanged (schema guaranteed upstream ⇒ haploid carries NA in hap2/both). Diploid output unchanged.
- **DISCOVERED (pass 5b, PENDING — blocks haploid validation):** `generate_summary_report.R` (run in `SUMMARY_REPORT`) has pervasive hap1/hap2 assumptions that will break/misrender on a haploid sample — `str_replace(id,"_hap[12]$","")` + `str_extract(id,"hap[12]$")` for sample/hap split (≥5 sites: 90, 164-165, 433-436), `across(c(hap1,hap2))` + `pivot_longer(c(hap1,hap2,both))` overview tables (191/200/270/488/495), `haplotype %in% c("hap1","hap2")` filters (241), and hardcoded `paste0(sid,"_hap1"/"_hap2")` + `for(hid in c(hap1_id,hap2_id))` loops in the contact-map / dotplot / riparian synteny section (348-371). Must be made ploidy-aware before the `ploidy=1` validation run reaches the report stage.

- **2026-07-03 — Phase 2, pass 3+4 (QUAST/MERQURY variable haplotype count).** `quast.nf` (QUAST) + `merqury.nf` (MERQURY) now take the ordered FASTA **list** + labels (QUAST) / meryl DB (MERQURY) instead of hardcoded `hap1_fasta, hap2_fasta` — one assembly (haploid) or two (diploid); QUAST copies its `QUAST_FINAL` list pattern, MERQURY builds `../asm` args dynamically. `assembly_qc.nf` re-pairing now emits `(sample_id, [fastas], [labels])` (labels `sample.<haplotype>`), feeds `QUAST(ch_paired_assemblies)` directly and `MERQURY` the `[fastas]`+db. Coupled — deploy together (module signatures changed). Diploid labels unchanged (`sample.hap1/.hap2`) ⇒ diploid QC identical; haploid now runs QUAST/MERQURY/BUSCO/MAPPING on one FASTA, but correct haploid **reporting** still needs pass 5 (R): the `primary` label lands in the `hap1` column there.

- **2026-07-03 — Phase 2, pass 1 (hifiasm ploidy branch + fork).** `hifiasm.nf`: `meta.n_hap`-driven — haploid (`n_hap==1`) forces `--primary` and disables Hi-C *phasing* (Hi-C still scaffolds downstream), emits one `<sample>.primary.p_ctg.fasta`; diploid unchanged. Output `assemblies` now variable-arity (glob `{hap1,hap2,primary}.p_ctg.fasta`); `gfa_hap1`/`gfa_hap2` (no consumers) consolidated into one globbed `gfa` output. `main.nf` fork (~566): destructure `(meta, fastas)`, coerce to list + `sort{it.name}`, `transpose`-zip against `forkHaplotypeMeta(meta)` → haploid 1 tuple, diploid 2. Diploid path unaffected by construction; haploid QC comes online after passes 3–4. D1(b)/D2/D3 chosen.

_(Add dated entries as work lands. Newest first.)_

- **2026-06-29 — Phase 1 REAL-RUN PARITY ✅ (CMat_203 identical to pre-refactor).** Real cluster run of `Sde-CMat_203` (diploid HiFi+Hi-C) through the refactored pipeline; QC metrics compared against the pre-refactor all-samples run. **All 153 compared values across 72 headline metrics matched exactly** (plus granular spot-checks: auN, error_rate, kmer_found identical to full precision), over all 8 stages (`ctg.base`→`final`), both haplotypes, and the combined `both` column. Includes the RNG-prone Hi-C mapping %/pairtools contact counts (identical → same inputs+resources) and the combined `both` QV (identical → `groupKey` regroup reproduces the old `groupTuple(size:2)` merge). QUAST total length + MERQURY kmer/QV identity ⇒ assembly FASTAs effectively confirmed unchanged. **Meta-map refactor is behavior-preserving on the diploid path.** Phase 1 done (threading + stub + parity). Remaining: §9 cleanup (non-blocking); Phases 2/4/5.
- **2026-06-29 — Phase 1 FULL LOCAL STUB PASS ✅ (Status: SUCCESS).** Every process in the meta-threaded DAG executed via stubs end-to-end, one sample (`Sde-CMat_203`, diploid). Execution-level confirmations: **per-hap fork** — all per-hap processes ran `2 of 2` (hap1+hap2) through mito-filter, misassembly correction, both FCS steps, Hi-C map/filter/metrics, scaffold r1+r2, gap-fill, teloclip, finalize, `CONTACT_MAP_FINAL`, compartments/TADs, TIDK, BUSCO, mapping-QC, snail, coverage; **per-sample re-pairing via `groupKey(sample,n_hap)`** — `QUAST`/`MERQURY`/`COMBINE_ASSEMBLY_QC` ran `1 of 1` with per-sample tags (`Sde-CMat_203:contig` … `:final`); **within-sample pairwise** — `PAIRWISE_ALIGNMENT (…_hap1_vs_…_hap2)`; **round-2 gating** fired; output shapes flowed to `SUMMARY_REPORT`. Validates wiring + tuple shapes + fork/regroup + gating end-to-end. NOT yet validated (needs the real run): tool output *content* and *parity* with `main`. Stub scaffolding (`stub_local.config` + the added stubs) kept for regression use.
- **2026-06-29 — Phase 1 stub fix: `COMPILE_FINAL_QC` stale stub.** Local stub reached the final compilation, then failed: `COMPILE_FINAL_QC`'s stub created `final_qc_report.tsv` (old, commented-out output) while the process declares `assembly_qc_metrics.csv` + `*.png` → "Missing output file" then real-`Rscript` 127. Pre-existing (original module, output block updated but stub not). Fixed the stub to touch `assembly_qc_metrics.csv` + `assembly_qc_metrics.png`. Audited the remaining terminal reporting stubs (`ASSEMBLY_REPORT`, `SNAIL_PLOT`, `COVERAGE_BOOK`, `SUMMARY_REPORT`) — all correct (touch their declared outputs). Expected to be the last stub fix for a full local `-stub` pass.
- **2026-06-29 — Phase 1 stub-coverage additions.** Local `-stub` execution (login node, `conda`/`singularity` disabled) confirmed the assembly spine + decontam DB setup stub through, then surfaced modules with **no `stub:` block** (pre-existing — absent in originals too; these always ran real commands under `-stub`). Added stubs (touch the declared outputs) to the process modules that run in the default config: `fcs_gx_screen.nf`, `fcs_clean_genome.nf`, `fcs_adaptor.nf`, `hic_mapping_metrics.nf` (both `HIC_BAM_METRICS` + `HIC_PAIRS_METRICS`), `hic_compartments.nf`, `hic_tads.nf`. Still stub-less but harmless by default: `fcs_dbgx.nf` (`FCS_DB_GET`, `storeDir`-cached → skips), the blobtools evidence chain (`GENERATE_DECONTAM_EVIDENCE`, gated off). Local stub bypass config: `process.executor='local'`, `process.scratch=false`, inflated `executor.cpus/memory`, `conda.enabled=false`, `singularity.enabled=false`.
- **2026-06-29 — Phase 1 STUB DAG CONSTRUCTION ✅ (wiring validated).** `nextflow run . -stub` (bare, one sample `Sde-CMat_203`) built the **full DAG** — every process materialized (spine, all `ASSEMBLY_QC_*` aliases, both scaffold rounds, final contact-map stage, pairwise/riparian/TIDK, `GENERATE_DECONTAM_EVIDENCE`, `COMPILE_FINAL_QC`→`SUMMARY_REPORT`) with correct per-sample tags, and SLURM scheduling began. No null-channel error, no `groupTuple` hang, no undefined-channel/shape failure. This validates the meta-map refactor at construction: channel definitions, tuple shapes, `join`/`groupKey`/`combine`/`mix` resolution, gating, per-hap fork. **Stub execution** (touch-through of all stubs to validate output-tuple shapes end-to-end) is blocked by an environmental issue, NOT the refactor: compute nodes lack conda init (`conda: command not found` in the task wrapper) — uniform across conda processes, while Singularity (`FCS_DB_GET`) and no-env (`FIND_MITO_REFERENCE`) tasks succeeded. Workaround for a bare stub: `-c stub.config` with `conda.enabled=false` + `singularity.enabled=false`; the real run gets conda via the launch wrapper.
- **2026-06-29 — Phase 1 stub bring-up fixes.** First `-stub` runs surfaced three non-threading issues, all fixed: (1) **`nextflow.config`** — added an `outdir` default to the `params{}` block; the reporting scopes referenced `params.outdir` before main.nf's line-38 default executed (config parses first). (2) **`functions/meta.nf`** — `buildMeta` now normalizes numeric ploidy (`1`/`2`) as well as `haploid`/`diploid`; the test samplesheet uses `ploidy=2`. (3) **`functions/parse_sample_sheet.nf`** — the per-row `catch` now logs the exception class + stack trace instead of a bare `${e.message}` (which printed `null`, masking an NPE thrown by a stale deployed `meta.nf`). Net: the stub now parses the sample sheet correctly. Remaining blocker is a data-path issue (Hi-C FASTQs absent at the samplesheet path), not pipeline code — parser is behaving correctly.
- **2026-06-29 — Phase 1, pass 7e (post-QC compilation — FINAL threading pass).** Rewired the reporting tail; all report modules unchanged (`COMPILE_FINAL_QC` / `ASSEMBLY_REPORT` / `SNAIL_PLOT` / `COVERAGE_BOOK` / `SUMMARY_REPORT`). Changes: (a) `COMPILE_FINAL_QC` metric maps `haplotype_id`→`meta` binding (cosmetic — only `tsv` extracted); (b) `SNAIL_PLOT_FINAL` input now `join`s `ch_finalized_assembly` × `ch_final_busco` on full `meta` (both per-hap), then maps to `meta.id` for the string-keyed module; (c) `COVERAGE_BOOK` input rebuilt — the old scalar `haplotype_id.replaceAll(/_hap[12]$/,'')` + `.combine(BAM_TO_FASTQ.out)` was doubly broken (map has no `replaceAll`; shape-mismatched combine now that `BAM_TO_FASTQ.out.fastq` is `(meta,fastq)`); replaced with the cross-granularity pattern (key on `meta.sample`, pair with per-sample HiFi, pass the **real** per-hap `meta` — `COVERAGE_BOOK` already accepted `meta`, the old site faked a minimal one); (d) manifest fixes — 5 entries read meta-keyed outputs (`FINALIZE_ASSEMBLY`, `CONTACT_MAP_FINAL`, `MITOHIFI`×2, `MITO_CIRCULAR_MAP`) and would have interpolated the whole map into the TSV `id` column — bound to `meta.id`/`meta.sample`. String-keyed viz outputs (SNAIL/PAIRWISE/RIPARIAN/TIDK) already emit clean ids, left as-is. **This completes Phase 1 — the workflow is meta-consistent end-to-end and ready for `-stub -resume`.** (This was the last `replaceAll` in `main.nf`.)
- **2026-06-29 — Phase 1, pass 7d (STEP 14 viz tail).** Added `ch_final_by_id` — a `(meta.id, fasta)` view of `ch_finalized_assembly` — and pointed the three leaf viz steps at it: pairwise/riparian (`toSortedList` source + both riparian `combine`s), `QUAST_FINAL` (label/assembly collection), and TIDK (`TIDK_EXPLORE`/`TIDK_SEARCH`). Modules unchanged; `PAIRWISE_ALIGNMENT`/`RIPARIAN_PLOT`/`QUAST_FINAL`/`TIDK_*`/`COLLECT_*` stay string-keyed. Pairwise `within_sample` `replaceAll(/_hap[12]$/,'')` retained deliberately — it runs on `meta.id`, correct for diploid pairing and (future) haploid no-pair. Contact-map / teloclip-QC / snail / coverage keep full `(meta, fasta)`. (Applied in two tries — the `QUAST_FINAL` `.map{list->…}.set{ch_quast_final_input}` tail was dropped on the first pass and restored.)
- **2026-06-29 — Phase 1, pass 7c (raw-read QC).** Threaded `hic_qc.nf` / `hifi_qc.nf` subworkflows and `fastqc_hic.nf` / `fastqc_hifi.nf` (per-hap `meta`, `meta.id == sample`); `multiqc_hic.nf` / `multiqc_hifi.nf` unchanged (aggregate-only). Rewired `HIC_QC_RAW` (reads-map extraction) and `HIFI_QC` (`.out.fastq`); `HIC_QC_TRIMMED` unchanged (already `(meta,r1,r2)`). Catalogued dead includes to §9.
- **2026-06-29 — Phase 1, pass 7b (QC call-sites + qc_mode).** Rewired all eleven `ASSEMBLY_QC_*` calls to pass their per-hap `(meta, fasta)` channel directly (`ASSEMBLY_QC_INITIAL` now sources `ch_contigs`, not `HIFIASM.out.assemblies`), deleting the eight `groupTuple(by:0,size:2)` prep blocks (~150 lines). Added `params.qc_mode` (`all_stages` default | `final_only`) via `run_all_qc`, gating the intermediate stages; the summary mix now builds from `Channel.empty()` respecting both `qc_mode` and which optional steps ran. No `groupTuple(by:0,size:2)` remains in `main.nf`.
- **2026-06-29 — Phase 1, pass 7a (assembly-QC core).** Rewrote the `ASSEMBLY_QC` subworkflow to take the per-hap `(meta, fasta)` stream and re-pair internally via `groupKey(meta.sample, meta.n_hap)` — one place, replacing both the `groupTuple(by:0,size:2)` re-pair and the `replaceAll(/_hap[12]$/,'')`. Threaded `build_meryl_db.nf`, `busco.nf`, `mapping_qc.nf` (per-hap `meta`); `quast.nf`/`merqury.nf`/`combine_assembly_qc.nf` unchanged (sample-string keyed, from the re-pair). Verified against the production original that a `.set{}` per-hap channel can feed both a processing step and a QC call (`ch_hifiasm_output` → Inspector + `ASSEMBLY_QC_PURGED` simultaneously), so DSL2 forks `.set{}` variables — call sites can pass per-hap channels directly.
- **2026-06-29 — Phase 1, pass 6b (final contact-map).** Threaded `contact_map.nf`, `hic_compartments.nf`, `hic_tads.nf` (per-hap `meta.id`) and rewired STEP 14's `make_final_contact_maps` block: reads combined to per-hap finalized assemblies on `meta.sample`; `MAP_HIC_TO_FINAL`/`FILTER_HIC_BAM_FINAL` + `final_raw_map`/`final_filtered` checkpoints + `CONTACT_MAP_FINAL`; compartments/TADs consume `CONTACT_MAP_FINAL.out.mcool` (forked). All Hi-C checkpoints (contig → scaffold-round2 → final) now closed.
- **2026-06-29 — Phase 1, pass 6 (gap-fill / teloclip / finalize).** Threaded `gap_filling.nf`, `teloclip.nf` (`TELOCLIP_EXTEND`; `COLLECT_TELOCLIP_STATS` left meta-agnostic — it only aggregates staged files), `finalize_assembly.nf` (per-hap `meta.id`; deterministic scaffold ordering preserved). Rewired STEP 13 gap-fill combine + STEP 13b teloclip `if/else` to `meta.sample` combines; `FINALIZE_ASSEMBLY` call unchanged. `ch_finalized_assembly` is now the per-hap `(meta, fasta)` canonical channel feeding STEP 14 + QC. Logged the `GAP_FILLING` report-heredoc cosmetic issue to §9.
- **2026-06-29 — Phase 1, pass 5b (Hi-C scaffold round2).** Wired STEP 8.5–12 (Inspector-on-scaffolds via `CORRECT_MISASSEMBLIES_SCAFFOLD`, `DECONTAMINATE_ASSEMBLY_SCAFFOLD`, then `MAP_HIC_TO_SCAFFOLD`/`FILTER_HIC_BAM_SCAFFOLD`/`SCAFFOLD_HIC_ROUND2`) with round-2 checkpoints `scaffold_round2_raw_map` / `scaffold_round2_filtered` / `scaffold_round2_space`. Same `meta.sample`-combine / same-task `meta`-join pattern as the contig stage. STEP 9 unchanged (pure pass-through into the meta-ready scaffold decontam). Note: the `if (params.run_scaffold_round2)` wrapper (+ its `else` `Channel.empty()`) had to be re-added after the initial paste dropped it.
- **2026-06-29 — Phase 1, pass 5a (Hi-C contig stage).** Threaded the four shared Hi-C modules (`map_hic_to_assembly`, `filter_hic_bam`, `scaffold_hic`, `hic_mapping_metrics`; per-hap `meta.id`, with the TSV `haplotype_id` column header retained as the R-parser contract). Wired STEP 6–8: Hi-C reads combined to per-hap assemblies on `meta.sample`, same-task checkpoint joins on `meta`; contig checkpoints `contig_raw_map` / `contig_filtered` / `scaffold_space`. `ch_assemblies_for_qc` kept as a separate two-consumer channel to preserve the baseline fork topology. The four modules' process names are unchanged, so all `*_CONTIG`/`*_SCAFFOLD`/`*_FINAL` include-aliases keep working.
- **2026-06-29 — Phase 1, pass 4 (decontam).** Threaded the `DECONTAMINATE_ASSEMBLY` subworkflow (key now `meta`; consume-twice topology preserved) and its three FCS modules (`fcs_gx_screen`, `fcs_clean_genome`, `fcs_adaptor`, per-hap `meta.id`). `ch_individual_haplotypes` collapses from the `(haplotype_id, sample_id, fasta)` triple to per-hap `(meta, fasta)`. Deferred: `GENERATE_DECONTAM_EVIDENCE` (blobtools side-branch) and `ASSEMBLY_QC_CONTIG_DECONTAM` (QC pass).
- **2026-06-29 — Phase 1, pass 3 (purge + Inspector).** Threaded `purge_dups.nf` and `correct_misassemblies.nf` (per-haplotype `meta.id`). Collapsed the purge and Inspector split/re-pair onto per-haplotype `ch_mito_filtered`; `ch_hifiasm_output` and `ch_assemblies_for_decontam` now flow as per-haplotype `(meta, fasta)` (else-branch splits removed). Minor: corrected the `correct_misassemblies` stub (`.fa`→`.fasta`, added missing `small_scale_error.bed`).
- **2026-06-29 — Phase 1, pass 2 (mito + caching slim).** Slimmed `meta` — removed `hic_rounds`/`scaffolders` (see §3.1 CACHING). Threaded `mitohifi.nf`, `filter_mito_contigs.nf` (collapsed the `_HAP1`/`_HAP2` aliases into a single per-haplotype `FILTER_MITO_CONTIGS`), and `mito_circular_map.nf`. Rewired `main.nf` STEP 4a: fan the per-sample mitogenome over `ch_contigs` via `combine(by:0)` → emit per-haplotype `ch_mito_filtered`. Purge / Inspector / `ASSEMBLY_QC_MITO_FILTERED` consumers of `ch_mito_filtered` remain old-shape pending later passes.
- **2026-06-29 — Phase 1, pass 1 (assembly spine).** Added `functions/meta.nf` (`buildMeta` + `forkHaplotypeMeta` — canonical meta contract). Rewrote `functions/parse_sample_sheet.nf` as wide/header-driven, emitting `tuple(meta, reads)`. Threaded `meta` through `bam_to_fastq.nf`, `trim_hic.nf`, `hifiasm.nf`. Wired `main.nf` spine: parse → BAM_TO_FASTQ → TRIM_HIC → scalar-keyed combine → HIFIASM → single haplotype fork (`ch_contigs`). Branch not yet runnable — everything from mito-filtering onward is still old-shape and rethreads in the next pass.

---

## 7. Decisions & Open Questions

### Phase 2 open decisions (haploid) — pending Jason

- **D1 — What is the haploid assembly?** For `n_hap==1`, which FASTA becomes `primary`? (a) Take hifiasm's first/primary contigs (`hap1.p_ctg`) from the existing 2-fasta output and discard the second — minimal, no per-sample hifiasm mode change; or (b) run hifiasm in `--primary` mode for haploid samples (per-sample mode driven by `meta.n_hap`) — "more correct" for a true haploid/collapsed assembly, but `hifiasm_primary` is currently a global param, so this means per-sample assembler config. *Leaning (a) to start; (b) as a follow-up if the collapsed primary is biologically needed.*
- **D2 — Haploid label:** keep `primary` (accurate for collapsed) as `meta.haplotype`/output naming, vs relabel to `hap1` (flows through existing hap1 machinery with least R churn). *Leaning keep `primary` for file/output naming; handle the mapping only in the QC layer (D3).*
- **D3 — Haploid QC-CSV representation:** the schema is `hap1,hap2,both`. For a haploid sample, put the single assembly's metrics in the **`hap1` column, `hap2`/`both` = NA** (keeps the schema uniform across mixed haploid+diploid runs, minimal downstream disruption), vs add a dedicated `primary` column (schema churn). *Leaning `hap1` column + NA.* Note: MERQURY on one assembly yields no combined QV, so `both`=NA falls out naturally.

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
- [ ] **Mito reference for short-read:** `FIND_MITO_REFERENCE` (NCBI lookup by species) feeds MitoHiFi and currently runs unconditionally at pipeline start. Review whether MitoFinder needs its reference supplied differently, and whether this prep step must be adapted or gated for short-read / MitoFinder runs (Phase 4).

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

---

## 9. Deferred Minor Fixes / Cleanup Pass

Small, non-blocking issues found mid-refactor. None block Phase 1 threading or the `-stub`
validation; batch them into a dedicated cleanup pass before the real parity run (or whenever
convenient). Resolved items move to the Change Log.

- [ ] **`GAP_FILLING` report heredoc** (`modules/gap_filling.nf`): the report is written with a
  quoted heredoc (`cat > … <<'EOF'`), so the `$(date)` / `$(seqkit stats …)` / `$(grep -c …)`
  command substitutions land in the report as **literal text** instead of evaluated values.
  Groovy-interpolated names (`${meta.id}`, `${scaffold_fasta}`) render fine. Cosmetic (report file
  only) — the gap-filled FASTA is unaffected, and the stub is unaffected. Fix: unquote the
  delimiter (`<<EOF`), but first audit the body for any `$` that must remain literal.
- [ ] **Dead includes — remove** (`main.nf`): included but never invoked after the refactor.
  Verified by scanning every `include` name against its call sites. Delete the `include` lines
  (and the backing files, once confirmed not imported by any subworkflow):
    - `HIC_QC_FROM_BAM_RAW`, `HIC_QC_FROM_BAM_FILTERED` → `workflows/hic_qc_from_bam.nf`
    - `HIC_SCAFFOLD_QC` → `workflows/hic_scaffold_qc.nf`
    - (already commented out, safe to delete: the `SCAN_TELOMERES; COLLECT_TELOMERE_RESULTS`
      include + `modules/scan_telomeres.nf` — superseded by TIDK)
  Also check whether `hic_qc_from_pairs.nf` / `hic_mapping_qc.nf` are orphaned files (not imported
  anywhere) and remove if so. Non-blocking: dead includes don't run in `-stub` or real runs.
- [x] **`params.outdir` config parse-order** (`nextflow.config`): the reporting scopes (`timeline`/`report`/`trace`/`dag`, lines ~485–500) reference `${params.outdir}`, but `outdir` was only defined in `main.nf` (line 38), which runs *after* the config is parsed → `Unknown config attribute params.outdir` on any bare `nextflow run` (masked previously because the launch wrapper passes `--outdir`). Pre-existing (also on `main`), surfaced by the first `-stub` run. Fix: add `outdir = './results'` to the config `params {}` block; the `main.nf` line-38 assignment becomes a redundant no-op.
- [ ] **`GENERATE_DECONTAM_EVIDENCE` (blobtools evidence branch) — moved to §10 as its own project.** Non-functional even before this refactor; always run with `--decon_make_blobtools_evidence false`. The line-49 `replaceAll`-on-`meta` crash is only the first of what is likely a fully un-threaded branch, so it's a resurrection project, not a cleanup one-liner. Tracked in §10.

---

## 10. Improvements / Feature Backlog (post-Phase-2)

- [ ] **Resurrect the blobtools decontam-evidence branch (`GENERATE_DECONTAM_EVIDENCE`).** This whole side-branch (coverage + DIAMOND taxonomy + BlobTools2 plots/report on cleaned assemblies) is non-functional and has been for a while — always run with `--decon_make_blobtools_evidence false` (its gate `params.decon_make_blobtools_evidence` defaults `true` at main.nf:267, which is a trap worth flipping to `false` when this is picked up). It was deferred in Phase 1 pass 4 and never threaded to the meta model: `workflows/generate_decontam_evidence.nf:49` calls `haplotype_id.replaceAll(/_hap[12]$/,'')` on a `meta` LinkedHashMap. Scope for the project: thread the `take:` channels (`decontaminated`/`contaminants`/`action_reports`/`taxonomy_reports`, now `(meta, …)`) using `meta.sample` for the HiFi-reads join and `meta.id` for labels/output naming, audit its 5 modules (`MAP_READS_MINIMAP2`, `DIAMOND_BLASTX`, `BLOBTOOLS_CREATE`, `BLOBTOOLS_VIEWPLOT`, `FCS_BLOB_EVIDENCE_REPORT`) for the same string-id assumptions, then validate the branch actually produces sensible blob plots (it's never been run to completion post-refactor). Two call sites: main.nf:1327/1343.

- [ ] **Shared Hi-C readset across multiple samples in one run.** Today `FASTQC_HIC` names outputs after the raw read basenames, so pointing two samples at the *same* Hi-C files makes `MULTIQC_HIC` (raw and trimmed) fail with an input-filename collision at staging. Doesn't affect normal runs (distinct samples → distinct read files); only bites deliberate readset sharing (e.g. benchmarking one individual diploid-vs-haploid). Interim workaround for the Phase 2 validation: symlink the Hi-C to per-sample-distinct names. Proper fix: stage `FASTQC_HIC` inputs under sample-prefixed names (or have MULTIQC stage collected zips into per-sample subdirs) so aggregated names are unique regardless of shared inputs. (HiFi side already safe — `BAM_TO_FASTQ` renames output to `${meta.sample}.fastq.gz`.)
