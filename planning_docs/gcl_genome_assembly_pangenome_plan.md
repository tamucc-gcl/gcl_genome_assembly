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
| Current phase | **Pangenome build-out complete (2026-08-08).** A, B, C Tier 1, D (per-haplotype + per-individual PCoA/NJ), E, F, G (frozen interface + docs), MultiQC, and H (progressive, opt-in) all built and validated end-to-end on the 4-assembly test; citations complete. **Remaining (not build-out):** the deferred C Tier 2 pair (inversion cross-check + orthogonal SV validation) and the actual **n=10 production run** (bump cactus memory). Deferred: BUSCO-on-graph. |
| Built so far | `PANGENOME` subworkflow: `CACTUS_PANGENOME` (+gaf emit), `PANGENOME_STATS`, `PANGENOME_REF_FASTA`, `PANGENOME_VARIANTS` (+bcftools stats), `PANGENOME_MANIFEST`, `PANGENOME_GROWTH`, `PANGENOME_PLOTS` (+R), `PANGENOME_2D_VIZ` (scattered), `PANGENOME_QC`, `PANGENOME_ODGI_STATS_MQC`, `MULTIQC_PANGENOME`, `PANGENOME_PCA_NJ`+`PANGENOME_POPSTRUCT` (odgi-similarity PCoA/NJ, hap+individual), `PANGENOME_PROGRESSIVE`(+plot, opt-in), `PANGENOME_REPORT` (+R, +main-report §8). Plus `COLLECT_NAME_MAPS` (teloclip remap). Caching fixed (deterministic cactus input order). |
| Roadmap ahead | n=10 production run (bump cactus memory) → C Tier 2 (inversion + SV validation, coupled) as capacity allows. |
| Test data | **10 genomes** of *Spratelloides delicatulus*. All modules validated at n=4 haplotypes; the full n=10 gives the real growth/openness numbers. |
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
  only; contiguity gate (`pangenome_max_chrom_scaffold_mult`) **off by default** since
  2026-08-18 — Minigraph-Cactus needs only the *reference* chromosome-scale and drops
  unplaceable contigs *per contig* as ambiguous, so a fragmented-but-sound assembly is a real
  individual; **comparability gate** drops a collapsed (`n_hap == 1`) assembly when the
  species also has phased assemblies (`pangenome_allow_collapsed = false`); reference always
  kept; require ≥ `pangenome_min_haplotypes`; kept/dropped set logged with a
  per-assembly reason. PanSN naming: reference + same-individual siblings get flat full-id
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
  (`null` = contiguity gate off), `pangenome_cactus_extra` (''), `pangenome_use_gpu` (false); `singularity { }` block;
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
- [x] **Graph-derived ordination + NJ trees, per haplotype AND per individual** — built 2026-08-08.
      `PANGENOME_PCA_NJ` (`odgi similarity -D '#' -p 2 -d` → per-haplotype distances, PanSN grouping
      **confirmed** on the real graph: 4 haplotype groups incl. reference) → `PANGENOME_POPSTRUCT`
      (`pangenome_popstruct.R`: `cmdscale` PCoA + `ape::nj`) → **four** PNGs in report §8: haplotype
      PCoA + NJ (4 units, real at n=4) and individual PCoA + NJ (haplotypes aggregated to diploid
      individual by stripping PanSN hap markers; 2 units at n=4 → placeholder, real at n=10). Distance
      = `jaccard.distance` column (odgi outputs `*.distance` columns, used directly). Gated
      `pangenome_popstruct`. odgi (`cactus_tools`) + `r-ape` on `pairwise_alignment` — no plink.

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
- [x] **Main report inclusion** — validated 2026-08-08. `generate_summary_report.R` appends the
      fragment (sentinel-guarded `readLines`), numbered **§8 Pangenome** with Methods → §9, TOC entry
      added, all conditional on `has_pangenome` so non-pangenome runs render unchanged. Rendered
      correctly in `assembly_report.md` (graph/variants/openness/QC tables + figures resolving).

### G. Packaging for reuse — `[DONE]` (2026-08-08)
- [x] Froze the `PANGENOME` interface — full reference in `pangenome_interface.md` (params, gating,
      envs, outputs, ops notes).
- [x] Documented every `pangenome_*` param; consolidated params block in `pangenome_params_block.groovy`
      with the two previously-implicit params (`pangenome_popstruct`, `pangenome_progressive`) now
      **declared explicitly** — no more reliance on Nextflow `null` defaulting.
- [x] All container labels catalogued (`cactus_pangenome`, `cactus_tools`, `pangenome_variants`,
      `panacus`, `pangenome_popstruct`, + reused `multiqc`/`pairwise_alignment`/`summarize_assembly`).
- [x] Default-on vs opt-in confirmed: all analyses default on except `pangenome_progressive` (opt-in);
      whole module gated by `run_pangenome` (opt-in). Report degrades gracefully via `NO_*` sentinels.

### H. Progressive (incremental-construction) growth — **opt-in** — `[DONE]` (built 2026-08-08)
- [x] **`PANGENOME_PROGRESSIVE`** (`pangenome_progressive = false`, default off) — incremental
      `minigraph -cxggs` over the reference-first assembly order, adding one at a time and recording
      graph size (nodes/edges/bp) after each → empirical growth table + curve
      (`PANGENOME_PROGRESSIVE_PLOT`), folded into report §8. Complements the analytic panacus growth (E).
      **Built on minigraph, not a full cactus rebuild:** the SV-graph construction is the tractable
      progressive backbone; a base-level cactus rebuild per step (~N× the full build, needs the big
      node) remains the heavier alternative if the intermediate base-level graphs are ever needed.
      Reuses `ch_cactus_in` + `cactus_tools`/`pairwise_alignment`; no new env.

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

- **2026-08-18** — **Comparability gate added** (`pangenome_allow_collapsed`, default
  `false`). Turning the contiguity gate off admitted three kinds of assembly, only one of
  which belongs in the graph. **(1)** A phased diploid from HiFi+Hi-C that missed chromosome
  scale is a genuine haplotype pair whose only defect is contiguity — which
  `cactus-graphmap-split` already handles per contig. It belongs in. **(2)** A HiFi-only
  *primary* is a collapsed diploid: a phase mosaic, not a haplotype. It has no interpretable
  ordination position, breaks panacus's one-unit-per-haploid-genome assumption, and (with
  `run_purge_dups = false`) carries retained haplotigs that FINALIZE preserves —
  harmonization *demotes* contained haplotigs to `unplaced_N` but never removes a sequence
  (the `NROWS == NSEQ` check guarantees it), so they reach cactus and read as duplication.
  It stays out. **(3)** A short-read primary was already excluded by `!meta.shortread` and
  is unaffected. The gate is decided at the GROUP level — collapsed is dropped only when
  a phased peer of the same species exists — because `main.nf:198` strips `meta.ploidy`
  onto a sample-keyed side-channel to keep ploidy tweaks off the task hash, so `n_hap` is the
  only ploidy signal surviving to `PANGENOME`.

  **Route 1 (align a non-member assembly to the finished graph) is dropped from the plan.**
  It only fixes contiguity, which is case (1)'s problem — and case (1) does not need it,
  because building works. Cases (2) and (3) have a *representation* problem that
  assembly-to-graph alignment does not touch. Both belong in the read-based downstream route,
  where every individual is genotyped against the graph by the same procedure regardless of
  how its assembly turned out.

- **2026-08-18** — **Contiguity gate off by default.** `pangenome_max_chrom_scaffold_mult` 3 → `null`; the scaffold-count filter is now
  opt-in and the kept/dropped set is logged. It had been silently excluding both haplotypes
  of one individual from the *Sde* graph. Minigraph-Cactus requires only the *reference* to
  be chromosome-scale: non-reference contigs are assigned to reference chromosome components
  by minigraph alignment, and unplaceable contigs are dropped *per contig* as ambiguous by
  `cactus-graphmap-split`, which is the right granularity. A `*_primary` assembly is never a
  redundant second copy of an individual already in the graph (`forkHaplotypeMeta` returns
  either `[_primary]` or `[_hap1, _hap2]`, never both), so the artifact the gate was written
  to catch cannot arrive by that route.
  **Caveat carried forward:** graph-derived Jaccard distances, and therefore the PCoA and NJ
  tree in §8, are sensitive to per-haplotype contiguity — a haplotype with more
  clipped/ambiguous sequence looks more distant for assembly reasons, not biological ones.
  Emit per-haplotype bp-retained-in-graph and check it against PCoA axis 1 before reading the
  ordination biologically. `midpoint_root` in `pangenome_popstruct.R` roots at the midpoint of
  the longest tip-to-tip path, so a single inflated tip can capture the root.

- **2026-08-08** — **Workstream G done (interface frozen + documented).** `pangenome_interface.md`
  (params table, gating logic, compute-env/label table, output tree, ops notes) + a consolidated
  `pangenome_params_block.groovy`. The two params the code used only via Nextflow's null defaulting
  (`pangenome_popstruct` on, `pangenome_progressive` off) are now declared explicitly. This closes the
  pangenome build-out — remaining items are the deferred C Tier 2 pair and the n=10 production run.

- **2026-08-08** — **Workstream H (progressive growth) built, opt-in** + **citations updated**.
  `PANGENOME_PROGRESSIVE` (`cactus_tools`: incremental `minigraph -cxggs`, reference-first, one
  assembly added at a time) → empirical growth table (k, sample_added, nodes/edges/bp) →
  `PANGENOME_PROGRESSIVE_PLOT` (`pairwise_alignment`/ggplot2) → curve, folded into report §8 as
  "Progressive growth (empirical)" (complements the analytic panacus growth). Gated
  `pangenome_progressive` (default **false**); minigraph-level (a full-cactus rebuild per step is
  prohibitive). Reuses `ch_cactus_in`; no new env. Citations: added Minigraph-Cactus, minigraph, vg,
  odgi, panacus, vcflib/vcfbub, BCFtools, MultiQC, ape to `generate_summary_report.R` (`refs` + `keys`,
  pangenome tools gated on `has_pangenome`; MultiQC always-on). Initially designed on
  plink-over-the-VCF, then corrected: the deconstruct VCF has only **2 diploid samples** (CBau's two
  haplotypes collapsed to one; reference omitted) → would plot 2 points, not the 4 haplotypes wanted.
  Rebuilt graph-native: `PANGENOME_PCA_NJ` (`cactus_tools`: `odgi similarity -D '#' -p 2 -d`, grouped
  by PanSN haplotype → 4 units incl. reference) → `PANGENOME_POPSTRUCT` (`pairwise_alignment` + `r-ape`:
  `pangenome_popstruct.R`, `cmdscale` PCoA + `ape::nj`) → PCoA scatter + NJ tree, report §8. Non-degenerate
  at n=4. Gated `pangenome_popstruct` (default on). No plink env. (§8 Pangenome, Methods→§9, TOC
  entry, figures resolving). Adjacent report fix surfaced by the test: the teloclip §6 table showed
  pre-harmonization `scaffold_N` names (teloclip runs before `HARMONIZE_SCAFFOLDS`). Fixed via a new
  `COLLECT_NAME_MAPS` module (merges the per-assembly `*.harmonized_name_map.tsv` already carried in
  `HARMONIZE_SCAFFOLDS.out.assemblies`) routed through `REPORTING`→`SUMMARY_REPORT`→
  `generate_summary_report.R`, which now left-joins the teloclip contig to the map (`{hap}_scaffold_N`
  → `chrN_1`/`unplaced_N`), sentinel-guarded so non-harmonized runs are unchanged. Separately logged
  teloclip over-aggressiveness (126 extensions vs 22 tidk telomeres) as future-projects §E — its own fix.
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
