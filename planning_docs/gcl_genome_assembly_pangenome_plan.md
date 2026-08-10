# gcl_genome_assembly — Pangenome Incorporation Plan

**Purpose.** Control document for adding a **general-purpose pangenome subworkflow** to
`gcl_genome_assembly`: build a Minigraph-Cactus graph from the finalized, harmonized
assemblies of a species, describe and QC it, produce a self-contained report section, and
expose a clean set of files/handles that *other* pipelines (e.g. a short-read
mapping/popgen pipeline) can consume. The subworkflow must be liftable into other
pipelines with minimal coupling.

**How to use this doc.** Update the Status Dashboard and per-workstream checkboxes as work
lands. Record completed work in the Change Log with a date. When a decision changes, edit
§3 and note it in the Change Log. Keep committed in the repo so any session can resume.

**Status legend:** `[TODO]` not started · `[WIP]` in progress · `[DONE]` complete · `[DEFERRED]` parked
Task checkboxes: `- [ ]` open · `- [x]` done.

**Cross-reference:** the target stats/plots/QC and downstream-file set follow the pangenome
reporting research report (tiered description + QC + progressive analysis, and the
downstream index/file requirements). Report §7 maps each item to that report's Tier 1/2/3.

---

## Status Dashboard

| Field | Value |
|---|---|
| Current phase | **Workstreams A + B complete & validated (2026-08-08).** Full `PANGENOME` run finished end-to-end on the 4-assembly test (2 individuals: CMat_203 + CBau_104, both haplotypes) and produced the complete category-1/2/4 output set: curated+flattened graph bundle, reference FASTA, variant catalog, and downstream manifest. Built on the corrected **15-chromosome** reference (`chr1_1…chr15_1`, no spurious chr16). Workstreams C–H not started. |
| Built so far | `PANGENOME` (gate + PanSN + species label, emits all 4 output categories); `CACTUS_PANGENOME` (build + curate + flatten + GPU toggle + `raw_vcf` handle); `PANGENOME_STATS` (vg/odgi, clip); **`PANGENOME_REF_FASTA`** (vg paths -F -S → reference.fa+.fai); **`PANGENOME_VARIANTS`** (vcfbub → SNP/indel/SV catalog + SV size spectrum); **`PANGENOME_MANIFEST`** (role→file TSV); config params + `pangenome_variants` label; main.nf wiring. Upstream: harmonization relative-size guardrail + batch-consensus check (produces the clean 15-chrom reference). |
| Roadmap ahead | Workstreams E → C+D → F → H (see §5, build order §8). Progressive rebuild (H) is LAST and opt-in. |
| Test data | **10 genomes** of *Spratelloides delicatulus* (clupeid fish, ~1 Gb, 15 chromosomes). Current validation used a 4-haplotype subset; the growth/core/Heaps analysis (E) becomes meaningful at the full n=10. |
| Last updated | 2026-08-08 |
| Production branch | `main` (S. delicatulus runs — DO NOT break) |

---

## 1. Goal & Scope

**Goal.** Turn the raw graph into a *described, QC'd, report-ready, and downstream-consumable*
pangenome resource, as a self-contained subworkflow reusable in other pipelines.

**In scope:**
- Graph construction (done) + curation/flattening (done).
- Descriptive statistics, quality diagnostics, visualizations, and openness/growth analysis.
- A self-contained **Markdown report section** emitted standalone and insertable into the main report.
- A clean **downstream manifest + file/handle set** for external pipelines.
- **Optional** progressive (n=1..N rebuild) mode.

**Explicitly out of scope (kept general-purpose):**
- **Short reads.** No read mapping, genotyping, or popgen here. Those live in a *separate*
  pipeline that consumes this one's outputs (GBZ, hapl, snarls, reference FASTA, VCF) via the manifest.
- **Read-mapping indexes** (`.dist` / `.min`). Deferred to downstream — how the graph is
  indexed depends on the consuming project. We expose GBZ + `.hapl` + snarls + reference FASTA + VCF; the mapping pipeline builds its own sample-specific indexes (haplotype sampling).
- Demographic modelling, PCA-on-population, SFS — downstream (this pipeline only produces the graph-level variant catalog + graph-derived PCA/tree on the *assembly* haplotypes).

---

## 2. Current State (as built, 2026-08-06)

- **`workflows/pangenome.nf` — `PANGENOME`.** Take: `(meta, fasta, fai)` finalized assemblies,
  `(taxid, ref_id)` harmonization reference, `(taxid, species)` resolved name. Gate: long-read
  only; exclude assemblies with chromosome-scale scaffold count > `pangenome_max_chrom_scaffold_mult` ×
  reference chromosome count (drops `*_hifi_primary`); reference always kept; require ≥
  `pangenome_min_haplotypes`. PanSN naming: reference + same-individual siblings get flat full-id
  names; other individuals group `<individual>.<hap>`. Species-derived publish/output label.
  Emit: `gbz`, `vcf`, `stats`, `versions`.
- **`modules/cactus_pangenome.nf` — `CACTUS_PANGENOME`.** Builds two-column seqfile; derives
  `--refContigs` from the reference FASTA's `chrN_p` names; runs `cactus-pangenome --vcf --haplo
  --gfa/gbz/viz/odgi/chrom-og full clip`. HOME/TMPDIR pointed at the scratch task dir (Toil fix).
  Curation: removes `chrom-subproblems/` + `chrom-alignments/` scratch; keeps everything else.
  Publish flattened via `saveAs` (strips `out/`). Clip-graph handles: `gbz/og/hapl/snarls/gfa/
  vcf/vcf_tbi/chrom_og/viz` + `all` catch-all. GPU toggle (`-gpu` image, `--gpu 1`, `--nv`).
- **`modules/pangenome_stats.nf` — `PANGENOME_STATS`.** `vg stats -lz` + `odgi stats -S` on the
  clip graph (`cactus_tools` CPU label).
- **Config:** `run_pangenome` (false), `pangenome_min_haplotypes` (2), `pangenome_max_chrom_scaffold_mult`
  (3), `pangenome_cactus_extra` (''), `pangenome_use_gpu` (false); `singularity { }` block;
  `cactus_pangenome` (toggled GPU) + `cactus_tools` (CPU) labels.
- **`main.nf`:** `ch_finalized_with_fai` join; `PANGENOME(finalized+fai, reference_id, species)`.
- **Harmonization:** `HARMONIZE_SPECIES` emits `${taxid}.reference_id.txt`; `HARMONIZE_SCAFFOLDS`
  surfaces `reference_id`.

---

## 3. Locked Decisions

1. **Curves without rebuild.** Growth / core / accessory / Heaps'-law curves are computed from
   the **single finished GFA** with `panacus` (permutation-averaged) — no repeated cactus builds.
2. **Progressive rebuild = opt-in, default off, built LAST** (`pangenome_progressive = false`).
   Only for when the actual intermediate graphs / per-step construction stats are wanted.
3. **Report built in the subworkflow.** `PANGENOME` produces its own report section so any
   embedding pipeline gets it for free.
4. **Report format = Markdown.** Emitted standalone (own rendered doc) *and* insertable into the
   main `generate_summary_report.R` report. **Reuse the existing R report stack** (no Quarto).
5. **Report always on** when a pangenome is built (no separate gate).
6. **Indexing (`.dist`/`.min`) deferred to downstream.** Expose GBZ + `.hapl` + snarls +
   reference FASTA + VCF; the mapping pipeline builds its own indexes.
7. **All QC default on** — including BUSCO-on-graph and orthogonal SV validation.
8. **Short reads excluded** — separate pipeline consumes this one's outputs via the manifest.

---

## 4. Target End-State — `PANGENOME` output contract

Four clean output categories (this *is* the generalization work — freeze this interface):

1. **Graph bundle** — GBZ, GFA, OG, snarls, `.hapl`, per-chromosome OG, decomposed VCF, +
   **reference-path FASTA** (new). Curated/flattened. *(Graph products; A.)*
2. **Downstream manifest** — `pangenome_manifest.tsv` (path · role · graph label) enumerating the
   files an external pipeline consumes, so mapping/popgen ingest programmatically. *(A — the reuse seam.)*
3. **Stats** — tidy machine-readable tables (stats JSON/TSV per graph) + all plots. *(B, C, E.)*
4. **Report fragment** — self-contained Markdown section + its figures + stats JSON, emitted
   whether or not a main report exists. *(F.)*

---

## 5. Workstreams

### A. Downstream-ready outputs & manifest — `[DONE]` *(the reuse seam)* — validated 2026-08-08
- [x] Surface named handles through `PANGENOME` `emit:` (`gbz`, `gfa`, `og`, `snarls`, `hapl`,
      `vcf`, `chrom_og`, `viz`, `ref_fasta`, `variants`, `sv_sizes`, `manifest`, `stats`).
- [x] **Reference-path FASTA** — `PANGENOME_REF_FASTA`: `vg paths -x gbz -F -S <ref>` → faidx →
      `Spratelloides_delicatulus.reference.fa`(+.fai). Produced on the run. Downstream surjection seam.
- [x] **`pangenome_manifest.tsv`** — `PANGENOME_MANIFEST`: role · file · label over gbz/gfa/og/
      snarls/hapl/variants-vcf(+tbi)/reference-fa(+fai). Produced on the run.
- [x] Decision-6 confirmed: `.dist`/`.min` NOT built here (left to the downstream mapping pipeline).

### B. Descriptive statistics — Report Tier 1 — `[WIP]` *(variant catalog done 2026-08-08)*
- [x] **Variant catalog** — `PANGENOME_VARIANTS`: `vcfbub --max-level 0 --max-ref-length` from the
      raw VCF → `bcftools norm -m -any` → length-classify SNP / indel (<`pangenome_sv_min_bp`) /
      SV → `variant_summary.tsv` (counts), `sv_sizes.tsv` (SV size spectrum + INS/DEL/COMPLEX),
      `variants.vcf.gz`(+tbi, the filtered catalog the manifest advertises). Produced on the run.
      *(INV/DUP not separately typed — deferred realignment refinement.)*
- [ ] Parse `vg stats` / `odgi stats` (produced: `vg_stats.txt`, `odgi_stats.txt`) into the report
      stats table. **→ handled in the report step (F).**
- [ ] **Graph-size-vs-genome sanity** — graph bp vs `ESTIMATE_GENOME_SIZE` / assembly total.
      **→ handled in the report step (F).**
- [ ] Non-reference/novel sequence + core/accessory partitioning — **computed in E (panacus).**

### C. Quality diagnostics — Report Tier 2 (**all default on**) — `[WIP]` *(Tier 1 built 2026-08-08)*

**Tier 1 — graph-intrinsic** (`PANGENOME_QC`, `cactus_tools`, gated `pangenome_qc`) — `[DONE]` validated 2026-08-08:
- [x] **Acyclicity** — `vg stats -A` → cyclic (expected with inversions — informational).
- [x] **Node degree / depth** — `odgi degree -S` (avg 2.73 / max 8) / `odgi depth`.
- [x] **Tangling / linearity** — `odgi stats -l -s` after `odgi sort -O` (mean_links_length 81.46,
      sum_path_node_dist 9.42).
- [x] **Re-alignment edit rate** — cactus `.gaf`: edit_rate 0.051 / graph_identity 0.949 / realigned
      4.27 Gb. *(Minigraph-GAF fidelity — coarse; a base-level metric would re-align to the GBZ.)*
- [x] Flags/parsing finalized: `vg stats -A` OK; awk `%d`→`%.0f` (INT_MAX fix); `all_paths`-row parse
      for linearity; odgi `-m` needs separate flags + `odgi sort -O` (compacted node IDs).

**Pangenome MultiQC** (`PANGENOME_ODGI_STATS_MQC` + `MULTIQC_PANGENOME`, gated `pangenome_multiqc`)
— `[DONE]` validated 2026-08-08: odgi module (whole + 15 per-chrom `*.og.stats.yaml`) + bcftools
module (`bcftools stats` on the catalog) → `<taxid>_pangenome_multiqc.html`. Reuses existing envs.

**Tier 2 — extrinsic** (need inputs from other subworkflows + new tool envs) — `[TODO]`:
- [DEFERRED] **BUSCO on the flattened graph** — skipped for now (may revisit). Flattened graph inlines
  alt alleles (corrupts gene models); reference-path BUSCO is identical to the reference assembly
  BUSCO already computed → redundant. The per-haplotype assembly BUSCOs remain the authoritative
  completeness figures. *(If revisited: `odgi flatten` → BUSCO as a caveated sanity-check; needs the
  `ch_busco_db` lineage-map threaded into `PANGENOME` + `params.busco_downloads`.)*
- [ ] **Spurious-inversion cross-check** — graph inversions × scaffold boundaries × harmonization
      composite/overlap flags → a **flag table** (not auto-fix). Cleaner, non-redundant signal.
      Inversion source = reverse-strand blocks in reference-vs-haplotype alignments (either wire the
      `FINAL_VIZ` all-vs-all PAFs in, or recompute self-contained from `ch_finalized`).
- [ ] **SV fidelity vs. orthogonal callers** — SVIM-asm / SyRI on the all-vs-all minimap2 from
      `FINAL_VIZ`; intersect with graph SVs (concordance table). Needs the all-vs-all PAFs wired in +
      a caller env. Most work. *(Overlaps the inversion check — both consume the pairwise alignments.)*

**NB** both remaining Tier 2 items are a coupled unit: they share the reference-vs-haplotype pairwise
alignments and both need cross-subworkflow wiring (or duplicated alignment). Best tackled together.

### D. Visualizations — `[WIP]` *(2D layout + R figures built 2026-08-08)*
- [x] 1D `odgi viz` (produced by cactus — have).
- [x] **2D layout** — `PANGENOME_2D_VIZ`: per-chromosome `odgi layout` + `odgi draw` on the clip
      per-chr graphs (best-effort; gated by `pangenome_2d_viz`, default on). Label `cactus_tools`.
- [x] **R report figures** — `r_scripts/pangenome_plots.R` + `PANGENOME_PLOTS` (label
      `pairwise_alignment`, ggplot2): growth/core curves with Heaps'-law fit **γ** and a
      rarefaction (Coleman) confidence band + machine-readable `growth_fit.tsv`; coverage-histogram
      core/accessory/private U-curve; **SV size histogram**; variant-class bar. All derived from the
      coverage histogram + variant catalog (no fragile histgrowth parse). Verified rarefaction math
      offline (core(N)=h(N); band ≈ hairline at bp scale, as expected).
- [ ] **Graph-derived PCA + NJ tree** from the decomposed VCF — `[DEFERRED]`: degenerate at the
      current n (2 samples / 4 haplotypes), needs a plink/ape env, and is arguably downstream-popgen
      territory. Build at the full n=10, or hand to the mapping pipeline.

### E. Openness / growth — Report Tier 3 (default via panacus) — `[WIP]` *(process built 2026-08-08)*
- [x] **`PANGENOME_GROWTH`** — panacus on the finished clip GFA, `--groupby-haplotype`:
      `histgrowth` (pangenome + core + soft-core curves, **exact expected = permutation-averaged**)
      and `hist` (coverage histogram → core/accessory/private partition). Data only. Env: `panacus`
      label (`bioconda::panacus`, no plotting deps). Awaiting run.
- [ ] Confirm the panacus CLI flags on the installed version (`--groupby-haplotype`, `-c/-l/-q/-t`)
      and the hist TSV column layout for `core_accessory.tsv`.
- [ ] **Plotting + Heaps fit → workstream D (R):** styled growth/core plot with a thin rarefaction
      band, coverage-histogram core/accessory bar, and the machine-readable γ + open/closed call,
      all from the histgrowth TSV. (`panacus-visualize` is deprecated/removed-soon, so not used;
      panacus `report` HTML is an optional standalone extra if wanted later.)

### F. Report fragment + main-report integration — `[WIP]` *(fragment built 2026-08-08)*
- [x] **`PANGENOME_REPORT`** — `r_scripts/pangenome_report.R` (label `summarize_assembly`): reads
      qc_metrics + growth_fit + variant_summary + odgi graph stats → **self-contained markdown
      section** (`pangenome_report.md`, same `c(md, ...)` style as the main report) + **stats JSON**
      (`pangenome_stats.json`). Gated `pangenome_report` (default on); robust to disabled sub-analyses
      via `NO_*` sentinels. Image paths written relative to outdir root so they resolve in the main report.
- [x] **Standalone emit** — the fragment + JSON publish to `pangenome/<taxid>/`.
- [ ] **Main report inclusion** — 4 edits to apply (delicate/large files): `generate_summary_report.R`
      (`--pangenome_report` arg + `readLines` append before final write), `SUMMARY_REPORT` (input +
      arg), `REPORTING` (take + pass-through), `main.nf` (`PANGENOME.out.report` → single path /
      `NO_PANGENOME` sentinel → `REPORTING`). Pipelines without a pangenome are unaffected.

### G. Packaging for reuse — `[TODO]`
- [ ] Freeze the `PANGENOME` interface (take: finalized+fai, reference_id, species; emit: the four categories).
- [ ] Document every `pangenome_*` param (§6).
- [ ] Add new container labels (§6) alongside `cactus_pangenome` / `cactus_tools`.
- [ ] Confirm default-on vs opt-in per analysis (§6 param defaults).

### H. Progressive rebuild — **opt-in, LAST** — `[DEFERRED / TODO-last]`
- [ ] **`PANGENOME_PROGRESSIVE`** (`pangenome_progressive = false`) — loop n=1..N over the first n
      genomes, reusing `CACTUS_PANGENOME` + `PANGENOME_STATS`; plot per-step construction stats.
      ~N× build cost; wants the H100 / 1.9 TB node. Only for the actual intermediate graphs.

---

## 6. Parameters & environments (to add)

**New params (all `pangenome_*`; QC toggles default true per Decision 7, exposed so downstream
users can disable):**
- `pangenome_progressive = false` — opt-in n=1..N rebuild.
- `pangenome_sv_min_bp = 50` — SV size threshold for the variant catalog.
- `pangenome_vcfbub_max_ref = 100000` — `vcfbub` bubble-size filter.
- `pangenome_growth_permutations = <n>` — panacus sampling-order permutations.
- `pangenome_busco = true` · `pangenome_sv_validation = true` · `pangenome_pca = true` — QC/vis
  toggles (default on; present for general-purpose reuse).
- `pangenome_busco_lineage` — reuse the pipeline's resolved BUSCO lineage.

**New process labels (containers/envs):**
- `panacus` — `bioconda::panacus` (growth/core/Heaps).
- `pangenome_variants` — `bioconda::bcftools bioconda::vcflib` (+ vcfbub) (variant catalog).
- `pangenome_sv_validate` — `bioconda::svim-asm bioconda::syri` (+ truvari/bcftools) (SV validation).
- BUSCO — reuse existing BUSCO label.
- R plots/report — reuse existing R report env (`generate_summary_report` stack) + any extra
  R packages for graph PCA/NJ (e.g. `ape`, `SNPRelate`/`adegenet`).
- `odgi` / `vg` ops — reuse `cactus_tools` (both bundled in the cactus image).

---

## 7. Report contents → research-report tier mapping

| Report item | Tier | Tool | Workstream |
|---|---|---|---|
| Graph vital-stats table (nodes/edges/paths/bp, per-chr) | 1 | `vg stats`, `odgi stats` | B (have) |
| Graph size vs genome size | 1 | arithmetic + genome size | B |
| Variant summary (SNP/indel/SV counts) | 1 | vcfbub + bcftools stats | B |
| SV size spectra + type breakdown | 1 | decomposed VCF + R | B/D |
| Non-reference/novel seq; core/accessory U-curve | 1/3 | panacus | E/D |
| Acyclicity; node degree/depth; tangling/linearity | 2 | vg/odgi | C |
| Re-alignment edit rate | 2 | cactus logs / `.gaf` | C |
| BUSCO completeness/duplication (graph) | 2 | BUSCO | C |
| SV fidelity vs SVIM-asm/SyRI | 2 | svim-asm/syri + intersect | C |
| Spurious-inversion flag table | 2 | graph SVs × scaffold/flags | C |
| 1D `odgi viz`; 2D `odgi draw` | — | odgi | D (1D have) |
| Growth / core / Heaps curves + exponent | 3 | panacus | E |
| Graph-derived PCA + NJ tree | — | VCF → R | D |
| Progressive per-step stats (opt-in) | 3 | cactus loop | H |

---

## 8. Build order

1. **A + B** — manifest, reference FASTA, emits, variant catalog. *(Makes the graph consumable &
   described fastest; unblocks the downstream mapping pipeline; fewest open decisions.)*
2. **E (panacus)** — growth/core/Heaps curves. *(Headline n=10 result.)*
3. **C (all QC, default) + D (all vis)** — the former "extras" (BUSCO-on-graph, orthogonal SV
   validation, PCA/NJ) are now default and fold in here, moved up per decision.
4. **F** — report fragment + main-report integration. *(Ties it together.)*
5. **H — progressive rebuild** — opt-in, LAST.

---

## 9. Open decisions

- None blocking. To confirm during build: exact `panacus` subcommand/flags for the coverage
  histogram + Heaps fit; whether 2D `odgi draw` is whole-genome, per-chr, or both (cost); report
  fragment mechanism (knitr child vs. include) once the main report engine is re-examined.

---

## 10. Change Log

- **2026-08-08** — **Workstream F fragment built** (`PANGENOME_REPORT` + `pangenome_report.R`, label
  `summarize_assembly`, gated `pangenome_report`). Reads qc_metrics + growth_fit + variant_summary +
  odgi graph stats → self-contained markdown section (`pangenome_report.md`) + stats JSON
  (`pangenome_stats.json`), robust to disabled sub-analyses via `NO_*` sentinels; image paths
  relative to outdir root. Wired into `PANGENOME` (base on the always-present variant summary,
  remainder-join the rest with sentinel fallback). Main-report inclusion is 4 documented edits to
  `generate_summary_report.R` / `SUMMARY_REPORT` / `REPORTING` / `main.nf` (sentinel + `readLines`
  append), left for the user to apply since those files are delicate/large. `PANGENOME_QC`
  works: `vg stats -A` → cyclic (expected), edit_rate 0.051 / graph_identity 0.949, and the odgi
  linearity/degree output is correct. Fixed: (1) `realigned_bp` hit INT_MAX — awk `printf %d` 32-bit
  overflow → `%.0f`; (2) `mean_links_length`/`sum_path_node_dist` were NA — values sit on the
  `all_paths` row, not the header → proper awk parse, and added `avg_node_degree` (2.73) /
  `max_node_degree` (8). (3) The odgi `.og.stats.yaml` files were **empty** — this odgi build rejects
  bundled short flags, so `odgi stats -m -sgdl` produced nothing; separated to `-m -s -g -d -l`.
  NB the edit_rate/identity is the **minigraph-GAF** re-alignment (coarse fidelity ≈ base divergence
  + small indels), not full base-level graph accuracy — 94.9% is reasonable for 2 divergent
  individuals + fragmented assemblies; a stricter metric would re-align inputs to the full GBZ (heavier).
  `PANGENOME_ODGI_STATS_MQC` (`cactus_tools`) emits `odgi stats -m -sgdl` YAMLs in the exact format
  MultiQC's odgi module ingests — one whole-genome + one per chromosome (`*.og.stats.yaml`) → a
  per-chromosome comparison of nodes/edges/paths/acyclicity/self-loops/composition/linearity.
  `bcftools stats` on the filtered catalog added to `PANGENOME_VARIANTS` (MultiQC bcftools module).
  `MULTIQC_PANGENOME` (`multiqc` label) aggregates both into `<taxid>_pangenome_multiqc.html`.
  Reuses existing envs (multiqc/cactus_tools/pangenome_variants) — no new labels. Images (1D viz /
  2D draw) can be folded in later via MultiQC custom-content; BUSCO-on-graph will slot in once wired.
  *(Also: the `gaf` emit on `CACTUS_PANGENOME` was added last turn but not presented — that was the
  `No such property: gaf` compile error; the corrected module is now provided.)*
  `pangenome_qc`, `cactus_tools`): acyclicity (`vg stats`), node degree/depth (`odgi`), linearity
  (`odgi stats -l -s` after `odgi sort -O`), and **re-alignment edit rate** from the cactus `.gaf`
  (the headline fidelity metric). Metrics are best-effort with NA fallbacks + a `qc_raw.txt` dump so
  the process never fails on one uncertain flag; flags/parsing to be finalized against the first run.
  Added a `gaf` emit to `CACTUS_PANGENOME` (output-only change — does not bust the cactus cache).
  Tier 2 (BUSCO-on-graph, orthogonal SV validation, inversion cross-check) needs cross-subworkflow
  wiring + new envs — outlined as the next sub-increment.
  `odgi sort -O` (optimize node IDs) resolved the "graph not optimized" index error and layout
  now runs. But the single task looped all 15 chromosomes serially (~25 min/chrom → ~6 h) and hit
  the 4 h SLURM limit (exit 140). Refactored `PANGENOME_2D_VIZ` to take ONE chromosome; the
  subworkflow scatters over the clip per-chrom `.og`s so the layouts run in parallel (~25 min each),
  with `errorStrategy 'ignore'` + a 4 h per-task limit so a slow/oversized chromosome can't kill the
  run. PNGs collected via `groupTuple`. (2D layout stays a viz-step concern; disable for very large
  production runs if it dominates runtime.) `groupTuple()` emitted assemblies
  in task-completion order, so the `val names` / `path fastas` inputs hashed unstably. The flatMap now
  sorts them deterministically (reference first, then by name) — restores caching and makes the graph
  build reproducible. NB this changes the input order once, so cactus rebuilds a single time on the next
  run, then `-resume` skips it. Also added `odgi sort -O` (optimize node IDs — required by
  `odgi layout`/`draw`'s index) before layout in `PANGENOME_2D_VIZ` (a viz-step concern, not in the
  build). *(First tried `-p Ygs`; odgi requires the graph to be optimized before that pipeline's
  index can build, so `-O` is the correct pre-step.)* panacus `histgrowth`/`hist` ran
  (`--groupby-haplotype` → N=4), variant catalog + R figures produced. Coverage histogram: private
  474.8 / accessory 224.7 / core 736.4 Mb (= 1.436 Gb graph); pangenome grows 999→1436 Mb; Heaps
  γ≈0.26. Fixed one R bug: panacus emits a **coverage-0 row** which broke R's 1-indexed histogram
  vector (shifted core/private) — `pangenome_plots.R` now drops coverage 0; awk `core_accessory.tsv`
  was already correct. **Pending:** `PANGENOME_2D_VIZ` produced no PNGs (odgi layout/draw failed under
  the best-effort guard — flag check needed); `panacus hist` format noted (coverage 0..N, col1=cov,
  col2=bp, with a non-`#` `panacus hist` line + `count\tbp` header). `PANGENOME_2D_VIZ` (per-chromosome
  `odgi layout`+`draw`, best-effort, gated by `pangenome_2d_viz`). `r_scripts/pangenome_plots.R` +
  `PANGENOME_PLOTS` (reuses `pairwise_alignment` R/ggplot2 env): growth/core curves + Heaps γ +
  rarefaction confidence band + `growth_fit.tsv`, coverage-histogram U-curve, SV size histogram,
  variant-class bar — all from the coverage histogram + variant catalog. Rarefaction math verified
  offline. PCA/NJ tree deferred (degenerate at current n; downstream-flavored). Awaiting rerun. panacus `histgrowth` +
  `hist` on the finished clip GFA, grouped by haplotype (PanSN). Emits growth/core curves,
  coverage histogram, and the core/accessory/private partition (data only; plotting + Heaps fit
  deferred to the R step in D — `panacus-visualize` is deprecated so unused). Wired into
  `PANGENOME` (gated by `pangenome_growth`, default on). New `panacus` label (`bioconda::panacus`)
  + `pangenome_growth_{count,coverage,quorum}` params. Awaiting run. (CMat_203 +
  CBau_104, both haplotypes). Confirmed: output flattening (`saveAs` strips `out/`) and scratch
  curation (chrom-subproblems / chrom-alignments removed); `PANGENOME_REF_FASTA`, `PANGENOME_VARIANTS`,
  `PANGENOME_MANIFEST` all produced their outputs; clip-only `PANGENOME_STATS`; full graph bundle
  retained (`.hal`/`.gaf`/`.paf`/`.sv.gfa`/`.raw.vcf`/`.stats.tgz`) per the general-purpose decision.
  Graph built on the corrected **15-chromosome** reference. Marked A `[DONE]`, B variant-catalog done.
- **2026-08-08** — Upstream fix that unblocked the clean reference: added a **relative-size guardrail**
  (`--min-chrom-frac`, drop chromosome-set members below a fraction of the median chromosome length)
  and a **batch-consensus check** (demote a reference chromosome with no chromosome-scale support from
  any other assembly) to `harmonize_names.py`. Removed the spurious 3.75 Mb "chr16" from the reference
  (a near-equal-size-cliff artifact); all four assemblies now resolve to 15 chromosomes.
  *(Also diagnosed CBau chr4/chr7 as genuine fragmentation at collapsed repeats — not a reference
  duplication and not a harmonizer bug; no action needed for the graph.)*
- **2026-08-06** — Plan created. Basic build validated (taxid 373251, n=2). Decisions §3 locked:
  panacus-not-rebuild for curves; progressive opt-in + last; report in subworkflow, Markdown, reuse
  R stack, always-on; indexing deferred downstream; all QC default-on; short reads out of scope.
