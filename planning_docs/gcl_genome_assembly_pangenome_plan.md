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
| Current phase | **Basic build working (2026-08-06).** `PANGENOME` subworkflow builds a Minigraph-Cactus graph end-to-end (validated: taxid 373251, 2 haplotypes, minigraph 271,543 events / 55 inversions, `vg`/`odgi` stats produced, outputs curated + flattened). Finalization workstreams A–G not started. |
| Built so far | `PANGENOME` subworkflow (gate + PanSN + species label); `CACTUS_PANGENOME` (build + curate + flatten publish + GPU toggle); `PANGENOME_STATS` (vg/odgi stats, CPU); config params + singularity block + `cactus_pangenome`/`cactus_tools` labels; main.nf wiring; harmonization `reference_id` emit. |
| Roadmap ahead | Workstreams A→G (see §5). Build order §8. Progressive rebuild is LAST and opt-in. |
| Test data | **10 genomes** of *Spratelloides delicatulus* (clupeid fish, ~1 Gb, 15 chromosomes). At n=10 the growth/core/Heaps analysis becomes meaningful. |
| Last updated | 2026-08-06 |
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

### A. Downstream-ready outputs & manifest — `[TODO]` *(the reuse seam)*
- [ ] Surface named handles through `PANGENOME` `emit:` (`gbz`, `hapl`, `snarls`, `gfa`, `og`,
      `vcf`, `chrom_og`, `viz`) — already emitted by `CACTUS_PANGENOME`; add to the subworkflow emit block.
- [ ] **Reference-path FASTA** — new step: `vg paths -x graph.gbz -F` (reference sample paths) →
      `samtools faidx`. Needed downstream for surject-to-BAM and reference-based callers. Graph product.
- [ ] **`pangenome_manifest.tsv`** — path · role · label for every downstream-relevant file.
- [ ] Decision-6 confirmed: do NOT build `.dist`/`.min` here.

### B. Descriptive statistics — Report Tier 1 — `[TODO]`
- [ ] Parse `vg stats` / `odgi stats` (have) into the stats table (whole-genome + per chrom-og).
- [ ] **Graph-size-vs-genome sanity** — compare graph bp to genome size (`ESTIMATE_GENOME_SIZE` /
      assembly total); flag over-collapse / over-inflation.
- [ ] **Variant catalog** — new process: `vcfbub` filter the decomposed VCF → classify SNP / indel
      (<50 bp) / SV (≥50 bp) relative to the reference path → counts, **SV size spectra**, **SV type
      breakdown** (ins/del/inv/dup). Env: `bcftools` + `vcflib`/`vcfbub`.
- [ ] Non-reference/novel sequence + core/accessory partitioning — computed once in E (panacus).

### C. Quality diagnostics — Report Tier 2 (**all default on**) — `[TODO]`
- [ ] **Acyclicity** — `vg stats -A` (reference path acyclic/unclipped).
- [ ] **Node-degree & node-depth** distributions — `odgi degree` / `odgi depth` (collapsed-repeat / tangle detection).
- [ ] **Tangling / linearity** — `odgi stats` (mean links length, sum-of-path-node-distances) after sort/layout.
- [ ] **Re-alignment edit rate** — parse from cactus logs / `.gaf` (graph faithfully represents inputs).
- [ ] **BUSCO on the flattened graph FASTA** — completeness / duplication. Reuse the pipeline's BUSCO env.
- [ ] **SV fidelity vs. orthogonal callers** — run SVIM-asm / SyRI on the all-vs-all minimap2
      alignments already produced in `FINAL_VIZ`; intersect with graph SVs (concordance table).
- [ ] **Spurious-inversion cross-check** — cross-reference graph inversions (55 seen at n=2) against
      scaffold boundaries + harmonization composite/overlap flags; surface as a **flag table** (not auto-fix).

### D. Visualizations — `[TODO]`
- [x] 1D `odgi viz` (produced by cactus — have).
- [ ] **2D layout** — `odgi layout` + `odgi draw` (whole-genome and/or per-chr).
- [ ] Plot scripts (R, `r_scripts/` staged-script convention): **SV size histograms**,
      **core/accessory U-curve**, **growth / Heaps curves**.
- [ ] **Graph-derived PCA + NJ tree** from the decomposed VCF (bubble presence/absence). Default-on,
      tolerant of small n (meaningful at n=10).

### E. Openness / growth — Report Tier 3 (default via panacus) — `[TODO]`
- [ ] **`PANGENOME_GROWTH`** — `panacus histgrowth` / `growth` on the final GFA → growth curve,
      core curve, coverage histogram (core/accessory/private), **Heaps'-law exponent + permutation
      band**, non-reference-sequence-by-haplotype-count. One cheap process; the headline n=10 result.
      Env: `panacus`.

### F. Report fragment + main-report integration — `[TODO]`
- [ ] **`PANGENOME_REPORT`** — R/RMarkdown render (reuse R report stack) taking the stats tables +
      all figures → **self-contained Markdown section** + **stats JSON**. Always on when a pangenome is built.
- [ ] **Standalone emit** — the fragment renders/exports on its own (own `.md` / rendered doc).
- [ ] **Main report inclusion** — `reporting.nf` / `generate_summary_report.R` detects the fragment
      (sentinel like the `NO_TELOCLIP` pattern) and includes it (knitr child) if present; pipelines
      without a pangenome are unaffected.

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

- **2026-08-06** — Plan created. Basic build validated (taxid 373251, n=2). Decisions §3 locked:
  panacus-not-rebuild for curves; progressive opt-in + last; report in subworkflow, Markdown, reuse
  R stack, always-on; indexing deferred downstream; all QC default-on; short reads out of scope.
