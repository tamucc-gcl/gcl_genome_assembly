# gcl_genome_assembly — Future Projects

Deferred / standalone efforts, split out of the active refactor plan on 2026-07-07. The active
round (through full parameter unification) is **complete**, and **§A (short-read organelle) is now
done** as well. Remaining: §B scaffolding/linked-reads (Hi-C rounds done; linked-read deferred —
no data), §C blobtools (untouched), §D docs (now timely).

Ordering is not fixed; these are scoped so any one can be picked up on its own.

---

## A. Short-read organelle tooling — **DONE (2026-07-29)** ✅

The 4a `ORGANELLE_ASSEMBLY` stub is replaced by a real `ORGANELLE` subworkflow
(`workflows/organelle.nf`) — one entry point branching on read type, generalized to N
organelle types, converging on the `FILTER_ORGANELLE` contract.

**As built:**
- **HiFi branch** → `MITOHIFI` (mito assembly + annotation) + `MITO_CIRCULAR_MAP`. MitoHiFi
  remains the sole HiFi organelle tool (mito only).
- **Short-read branch** → `SHORTREAD_ORGANELLE` (GetOrganelle, `getorganelle.nf`) +
  `ORGANELLE_ANNOTATION` (MITOS2 → GenBank + circular/gene maps). GetOrganelle runs once per
  resolved organelle target (`-F ${ORG}`), emitting per-organelle `(meta, organelle_type, file)`
  with `{circular, linear, failed}` status.
- **Kingdom-driven target selection** — `organelleTypesFor()` (plant → plastid + mito; else
  mito), built into `ch_organelle_by_taxid` in `main.nf`; per-organelle GetOrganelle tuning via
  `getorganelleRecursionFor` / `getorganelleKmersFor` / `getorganelleCoverageFor`.
- **Nuclear stripping** — `FILTER_ORGANELLE` baits per sample on ALL organelle assemblies
  (mito + plastid) and removes those contigs from the nuclear assembly.
- Emits `assemblies` (→ FILTER_ORGANELLE bait), `annotation`, `stats`, `circular`, `gene_map`,
  `notes` (plant organelles), `versions`.

- [x] `ORGANELLE` real subworkflow (HiFi→MitoHiFi; SR→GetOrganelle; multi-organelle; filtering).
- [x] Kingdom-driven target selection (plant → mito+plastid; animal → mito) via `organelleTypesFor`.
      *(per-run override for ambiguous cases — add if a case arises.)*
- [x] Validated on the three-sample test: plant SR → mito + plastid; HiFi → MitoHiFi unchanged;
      nuclear assembly stripped of organelle contigs.
- [x] `generate_summary_report.R` §5 ("Mitochondrial Genome") — **kept mito-scoped by design; no
      generalization needed.** §5 reports only organelle genomes with reliable command-line
      annotation + circularization: the HiFi/MitoHiFi mitogenome and (via `ORGANELLE_ANNOTATION` /
      MITOS2) the short-read **animal** mito. Plant organelles (plastid + mito) are assembled by
      GetOrganelle and used for nuclear-contig filtering, but are **intentionally excluded from the
      report** — no CLI tool reliably guarantees their circularization/annotation, so surfacing
      partial/unvalidated stats would mislead. They ride the subworkflow's `notes` emission
      (assembled + filtered, not annotated/reported). The mechanism already enforces this: plant
      organelles never produce `mito_stats` rows, so §5 filters them out automatically.

---

## B. Scaffolding chain + linked reads (was "Phase 5") — *later project*

Goal: generalize the scaffolding stage into a per-round, evidence-routed chain and add linked-read support. Completes the input matrix (short-read+Hi-C, HiFi/short-read+TellSeq).

- [ ] Refactor Hi-C scaffolding into a `SCAFFOLD_HIC_ROUND` subworkflow; express round1/round2 as chained calls; add `params.hic_scaffold_rounds` (global) so N rounds are data-driven rather than hardcoded.
- [ ] `.branch{}` routing by evidence: linked → hic → none; `.mix()` to converge on `ch_final_assembly`.
- [ ] Finishing gates by data type: gap-fill (TGSGapCloser long-read / ABySS-Sealer short-read), teloclip, contact maps.
- [ ] **Linked-read scaffolding module** (ARCS/LINKS lineage, optionally Tigmint) emitting `(meta, fasta)`. **Build to contract; live validation DEFERRED — no linked-read data on hand.**
- [ ] Revisit redundans as a scaffolder / gap-closer option alongside linked-read + Hi-C; decide how to factor it in.
- [ ] Validation: HiFi+Hi-C 3-round run works; short-read+Hi-C scaffolds and gap-closes; HiFi+TellSeq / short-read+TellSeq deferred (contract-only until data exists).

Notes carried from the refactor plan: the meta model already carries `meta.hic` / `meta.hifi` for evidence-gating; short-read+Hi-C and short-read+TellSeq combos were explicitly deferred to this project because they need the scaffolding chain.

---

## C. Blobtools decontamination-evidence resurrection (`GENERATE_DECONTAM_EVIDENCE`) — *independent*

This whole side-branch (coverage + DIAMOND taxonomy + BlobTools2 plots/report on cleaned assemblies) is non-functional and has been for a while — always run with `--decon_make_blobtools_evidence false`. Its gate `params.decon_make_blobtools_evidence` defaults **`true`** at `main.nf:267`, which is a trap worth flipping to `false` when this is picked up.

Scope:
- [ ] Thread the `take:` channels (`decontaminated` / `contaminants` / `action_reports` / `taxonomy_reports`, now `(meta, …)`) to the meta model — `meta.sample` for the HiFi-reads join, `meta.id` for labels/output naming. `workflows/generate_decontam_evidence.nf:49` still calls `haplotype_id.replaceAll(/_hap[12]$/,'')` on a `meta` LinkedHashMap (the crash).
- [ ] Audit its 5 modules (`MAP_READS_MINIMAP2`, `DIAMOND_BLASTX`, `BLOBTOOLS_CREATE`, `BLOBTOOLS_VIEWPLOT`, `FCS_BLOB_EVIDENCE_REPORT`) for the same string-id assumptions.
- [ ] Validate the branch actually produces sensible blob plots — it has never been run to completion post-refactor. Two call sites: `main.nf:1327/1343`.

Independent of everything else; do whenever blob plots are wanted.

---

## D. README / documentation revamp — *independent, after the shape settles*

A dedicated effort to fully revamp the README and write end-user + developer docs, once the input matrix and architecture stabilize (i.e. after the scaffolding/linked-read work at least). Not refactor cleanup — its own project.

### Scope
- [ ] Rewrite README: what the pipeline does, quick-start, install/deps (conda vs singularity, Nextflow 23.10.1), how to run on Crest/SLURM, the launcher wrapper.
- [ ] Input docs: the wide header-driven sample sheet — every column (required vs optional), worked examples for each input combo (HiFi+Hi-C diploid, HiFi haploid, HiFi-only, short-read-only, +TellSeq later).
- [ ] Parameter reference: full list grouped by stage, defaults, meaning (pairs with the unified param scheme delivered in 4b-i).
- [ ] Outputs/results layout: directory tree, the QC report, per-stage artifacts, the assembly webshare.
- [ ] Architecture / developer docs: the meta-map model, subworkflow map, how to add an assembler/scaffolder/organelle tool, **caching rules (what re-hashes and why — incl. the script-staging pattern from Track 1)**.
- [ ] Per-tool notes + citations.

## E. teloclip over-aggressive extension — *independent, small*

*Noted 2026-08-08 from the S. delicatulus pangenome-integration test report. Separate from the
teloclip **naming** issue (scaffold_N vs harmonized chrN_1), which is fixed via the harmonization
name-map remap in `generate_summary_report.R`.*

teloclip extends a scaffold end on **any** soft-clip containing a telomeric motif, so it fires on
small repeat-rich unplaced scaffolds that are not real telomeres. On the test (Sde-CBau_104 hap1):

- teloclip: **126 contigs extended / 172 extensions / 1.71 Mb added** (mean 9.9 kb, max 28.1 kb).
- tidk on the final assembly: only **22 scaffolds with a detectable telomere array** (12.79%).

So ~100 ends were extended on weak/spurious signal, adding ~1.7 Mb of read overhang, much of it
non-telomeric. §6 of the report therefore reads as contradictory: the teloclip subsection ("126
extended") vs the tidk subsection ("22 with telomeres"). This is teloclip's default permissiveness,
**not** a harmonization or pipeline bug.

### Options to scope (pick per appetite)
- [ ] **Report framing (cheapest — do regardless):** relabel §6 so "extension attempts" (teloclip)
      and "confirmed telomere arrays" (tidk) are clearly distinct metrics, not read as one number.
- [ ] **Tighten teloclip:** require a minimum telomeric-repeat run length / motif count in the
      soft-clip before extending, and/or restrict extension to chromosome-scale scaffolds (length
      filter) rather than all unplaced fragments.
- [ ] **Cross-gate with tidk:** accept a teloclip extension only if tidk confirms a telomere array
      at that end after extension — turns teloclip into propose-then-validate.
- [ ] **Retention decision:** decide whether the ~1.7 Mb of unconfirmed extension sequence should be
      trimmed (affects final length ~0.15%) or kept.

Scoped small; independent of the pangenome work and of the naming fix.

### Running notes to fold into the docs (from the refactor plan)
- **Parameter scheme + column↔param duality**: the unified grouped params (from 4b-i), and that optional strategy columns (`ploidy`, `n_hap`, `assembler`, `dedup`, `mito_tool`, organism identity) can be set globally as params or per-row, per-row winning.
- **`ploidy` (organism) vs `n_hap` (output haplotype count)**: distinct; the `n_hap` override; e.g. collapsed diploid = `ploidy=2` + `n_hap=1` (hifiasm `--primary`), NOT a false haploid.
- **Assembler selection** (`meta.assembler`: hifiasm | spades) and how short-read routes through the single collapsed/`primary` path.
- **Genome-size estimation** (jellyfish → GenomeScope2): `-p` = organism ploidy; where the estimate lands in the report.
- **Organelle handling**: MitoHiFi on HiFi (mito only), GetOrganelle on short-read; plant → mito + chloroplast via taxid → kingdom.
- **Evidence-gating**: which steps run for which input types (Hi-C chain needs `meta.hic`; long-read finishing needs `meta.hifi`).
- **Caching / script-staging**: scripts are declared `input: path(...)` so edits re-hash the dependent process (project `bin/` contents are NOT hashed — the gotcha that motivated Track 1).
- (append more as projects land)
