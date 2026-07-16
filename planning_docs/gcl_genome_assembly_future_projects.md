# gcl_genome_assembly — Future Projects

Deferred / standalone efforts, split out of the active refactor plan (`gcl_genome_assembly_refactor_plan.md`) on 2026-07-07. The **active** round finishes at Phase 4a (done) → Track 1 (script cache-staging) → **4b-i (identity/taxonomy consolidation + full parameter unification)**, which is the intended stopping point. Everything below is **after** that round — some soon (organelle tooling), some genuinely later (scaffolding/linked reads), some independent (blobtools, docs).

Ordering is not fixed; these are scoped so any one can be picked up on its own.

---

## A. Short-read organelle tooling (fills the 4a `ORGANELLE_ASSEMBLY` stub) — *next-ish*

Goal: replace the non-HiFi organelle no-op stub with a real subworkflow. **Decided:** MitoHiFi stays the **sole** organelle tool on the **HiFi branch** (mito only, until a HiFi chloroplast tool is identified); **GetOrganelle** (https://github.com/kinggerm/getorganelle) handles the **short-read-only** branch — **mito** for animals, **mito + chloroplast** for plants.

Depends on 4b-i (needs taxid → kingdom to decide mito-only vs mito+chloroplast, and taxid → name for references).

- [ ] `ORGANELLE_ASSEMBLY` (real): HiFi branch → MitoHiFi (as today); SR branch → GetOrganelle (assembly → annotation → circular map → organelle-contig **filtering** from the nuclear assembly, generalized to N organelle types). Both converge on the same downstream contract the 4a stub defined.
- [ ] Kingdom-driven target selection (plant → mito + chloroplast; animal → mito), from the taxid→kingdom lookup built in 4b-i. Per-run override for ambiguous cases.
- [ ] Extend `generate_summary_report.R` §5 ("Mitochondrial Genome") to report multiple organelles (mito + chloroplast) per sample.
- [ ] Validation: plant SR run yields mito + chloroplast; fish SR run yields mito; HiFi run unchanged (MitoHiFi); nuclear assembly correctly stripped of all organelle contigs.

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

### Running notes to fold into the docs (from the refactor plan)
- **Parameter scheme + column↔param duality**: the unified grouped params (from 4b-i), and that optional strategy columns (`ploidy`, `n_hap`, `assembler`, `dedup`, `mito_tool`, organism identity) can be set globally as params or per-row, per-row winning.
- **`ploidy` (organism) vs `n_hap` (output haplotype count)**: distinct; the `n_hap` override; e.g. collapsed diploid = `ploidy=2` + `n_hap=1` (hifiasm `--primary`), NOT a false haploid.
- **Assembler selection** (`meta.assembler`: hifiasm | spades) and how short-read routes through the single collapsed/`primary` path.
- **Genome-size estimation** (jellyfish → GenomeScope2): `-p` = organism ploidy; where the estimate lands in the report.
- **Organelle handling**: MitoHiFi on HiFi (mito only), GetOrganelle on short-read; plant → mito + chloroplast via taxid → kingdom.
- **Evidence-gating**: which steps run for which input types (Hi-C chain needs `meta.hic`; long-read finishing needs `meta.hifi`).
- **Caching / script-staging**: scripts are declared `input: path(...)` so edits re-hash the dependent process (project `bin/` contents are NOT hashed — the gotcha that motivated Track 1).
- (append more as projects land)
