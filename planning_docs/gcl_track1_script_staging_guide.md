# Track 1 — Stage the QC/report/plot scripts as process inputs (cache-correct)

**Goal:** make edits to the R/Python scripts actually invalidate the cache, and move them into `bin/r_scripts/` + `bin/python/`.

**Mechanism (the important bit):** Nextflow's task hash = session ID + process name + **command string** + **declared input files** (path/timestamp/size) + container/env + directives. Today the scripts are called as `Rscript ${projectDir}/r_scripts/compile_qc.R …`, which bakes an absolute path into the command string but leaves the script file **undeclared** — so editing it changes nothing in the hash and `-resume` skips the process. Putting the file under `bin/` does **not** fix this (project `bin/` is on `PATH` and bundled, but its *contents* aren't hashed). The fix is to declare each script as a process `input: path(...)` and call it via the staged variable — **staged inputs are hashed**, so an edit re-runs the process.

Because we stage each script explicitly, we do **not** rely on `bin/` being on `PATH`. The folder is purely organizational — `bin/r_scripts` + `bin/python` is fine (so would be `scripts/`).

---

## Step 1 — Move the scripts

```bash
mkdir -p bin/r_scripts bin/python
git mv r_scripts/*.R   bin/r_scripts/     # 5 files   (plain `mv` if not under git)
git mv py_scripts/*.py bin/python/        # 5 files
rmdir r_scripts py_scripts 2>/dev/null || true
```

- **→ `bin/r_scripts/`** (5): `compile_qc.R`, `combine_individual_assembly_qc.R`, `generate_summary_report.R`, `dotplot_paf.R`, `riparian_paf.R`
- **→ `bin/python/`** (5): `generate_assembly_report.py`, `bigwig_genome_book.py`, `plot_mito_circular.py`, `make_tad_book.py`, `plot_compartments_pc1_genomewide.py`

(`chmod +x` is optional — we invoke them via `Rscript`/`python3`, not by PATH lookup.)

---

## Step 2 — Declare the script files once (as `file()` value channels)

In `main.nf`, just after the `include { … }` block (~line 525), add:

```groovy
// ── Track 1: helper scripts declared as inputs so edits invalidate the cache ──
ch_compile_qc_script      = file("${projectDir}/bin/r_scripts/compile_qc.R",                     checkIfExists: true)
ch_summary_report_script  = file("${projectDir}/bin/r_scripts/generate_summary_report.R",        checkIfExists: true)
ch_dotplot_script         = file("${projectDir}/bin/r_scripts/dotplot_paf.R",                    checkIfExists: true)
ch_riparian_script        = file("${projectDir}/bin/r_scripts/riparian_paf.R",                   checkIfExists: true)
ch_assembly_report_script = file("${projectDir}/bin/python/generate_assembly_report.py",         checkIfExists: true)
ch_coverage_book_script   = file("${projectDir}/bin/python/bigwig_genome_book.py",               checkIfExists: true)
ch_tad_book_script        = file("${projectDir}/bin/python/make_tad_book.py",                    checkIfExists: true)
ch_compartments_script    = file("${projectDir}/bin/python/plot_compartments_pc1_genomewide.py", checkIfExists: true)
```

`checkIfExists: true` fails the run immediately if a script is missing/misnamed after the move — a cheap guard against typos.

The **other two** scripts are used inside subworkflows, so declare them *there* (shown in the table):
- `combine_individual_assembly_qc.R` → inside `assembly_qc.nf` (used by `COMBINE_ASSEMBLY_QC`)
- `plot_mito_circular.py` → inside `organelle_assembly.nf` (used by `MITO_CIRCULAR_MAP`)

---

## Step 3 — The per-process edit (uniform, three parts)

For every process: **(a)** add `path(<var>)` as the **last** line of its `input:` block; **(b)** change the script-invocation line to the staged variable; **(c)** add the `file()` ref as the **last** call-site argument. Keeping the script last everywhere means positional arg-matching stays trivial.

### Worked example A — `COMPILE_FINAL_QC` (plain `path` inputs) — `modules/compile_final_qc.nf`

```groovy
    input:
    path(assembly_summaries)
    path(bam_metrics)
    path(pairs_metrics)
    path(compile_qc_script)          // ← ADD
```
```diff
-    Rscript ${projectDir}/r_scripts/compile_qc.R \
+    Rscript ${compile_qc_script} \
```
Call site — `main.nf` ~1634:
```groovy
    COMPILE_FINAL_QC(
        ch_all_assembly_summaries.map { sample_id, qc_label, tsv -> tsv }.collect(),
        ch_all_bam_metrics.map { meta, checkpoint, tsv -> tsv }.collect(),
        ch_all_pairs_metrics.map { meta, checkpoint, tsv -> tsv }.collect(),
        ch_compile_qc_script          // ← ADD
    )
```

### Worked example B — `COVERAGE_BOOK` (tuple input) — `modules/coverage_book.nf`

```groovy
    input:
    tuple val(meta), path(assembly), path(hifi_reads)
    path(coverage_book_script)       // ← ADD
```
```diff
-    python3 ${projectDir}/py_scripts/bigwig_genome_book.py \
+    python3 ${coverage_book_script} \
```
Call site — `main.nf` ~1676:
```diff
-    COVERAGE_BOOK(ch_coverage_book_input)
+    COVERAGE_BOOK(ch_coverage_book_input, ch_coverage_book_script)
```
> A standalone `path(script)` input **alongside** a `tuple` input is a value channel → it broadcasts to every sample. That's exactly what we want (same script for all).

`${compile_qc_script}` etc. render as the **staged basename** (`compile_qc.R`) in the command — no `${projectDir}` — because staged inputs land in the task's working dir.

---

## All ten processes — apply the same three edits

| Process | Module | Invocation line: change | Call site → add last arg | `file()` decl in |
|---|---|---|---|---|
| `COMPILE_FINAL_QC` | `compile_final_qc.nf` | `Rscript ${projectDir}/r_scripts/compile_qc.R` → `Rscript ${compile_qc_script}` | `main.nf`~1634 (4th arg) `ch_compile_qc_script` | main.nf |
| `COVERAGE_BOOK` | `coverage_book.nf` | `python3 …/bigwig_genome_book.py` → `python3 ${coverage_book_script}` | `main.nf`~1676 `COVERAGE_BOOK(ch_coverage_book_input, ch_coverage_book_script)` | main.nf |
| `SUMMARY_REPORT` | `summary_report.nf` | `Rscript …/generate_summary_report.R` → `Rscript ${summary_report_script}` | `main.nf`~1819 (7th arg) `ch_summary_report_script` | main.nf |
| `ASSEMBLY_REPORT` | `assemblyReport.nf` | `python3 …/generate_assembly_report.py` → `python3 ${assembly_report_script}` | `main.nf`~1641 `ASSEMBLY_REPORT(COMPILE_FINAL_QC.out.metrics, ch_assembly_report_script)` | main.nf |
| `PAIRWISE_ALIGNMENT` | `pairwise_alignment.nf` | **line ~156 only** `Rscript …/dotplot_paf.R` → `Rscript ${dotplot_script}` (leave the inline `Rscript -e` at ~41) | `main.nf`~1242 `PAIRWISE_ALIGNMENT(ch_pairwise_input, SETUP_PAFR.out.ready, ch_dotplot_script)` | main.nf |
| `RIPARIAN_PLOT` | `riparian_plot.nf` | `Rscript …/riparian_paf.R` → `Rscript ${riparian_script}` | `main.nf`~1256 `RIPARIAN_PLOT(ch_riparian_input, ch_riparian_script)` | main.nf |
| `HIC_COMPARTMENTS` | `hic_compartments.nf` | `python …/plot_compartments_pc1_genomewide.py` → `python ${compartments_script}` | `main.nf`~1173 (add `ch_compartments_script` after `…max_contigs ?: 30`) | main.nf |
| `HIC_TADS` | `hic_tads.nf` | `python …/make_tad_book.py` → `python ${tad_book_script}` | `main.nf`~1180 (add `ch_tad_book_script` after `…max_contigs ?: 0`) | main.nf |
| `COMBINE_ASSEMBLY_QC` | `combine_assembly_qc.nf` | `Rscript …/combine_individual_assembly_qc.R` → `Rscript ${combine_qc_script}` | `assembly_qc.nf`~144 `COMBINE_ASSEMBLY_QC(ch_all_qc_labeled, combine_qc_script)` | **assembly_qc.nf** † |
| `MITO_CIRCULAR_MAP` | `mito_circular_map.nf` | `python3 …/plot_mito_circular.py` → `python3 ${mito_circular_script}` | `organelle_assembly.nf`~57 `MITO_CIRCULAR_MAP(MITOHIFI.out.annotation, mito_circular_script)` | **organelle_assembly.nf** † |

† In each subworkflow's `main:` block, before the call, add:
```groovy
combine_qc_script   = file("${projectDir}/bin/r_scripts/combine_individual_assembly_qc.R", checkIfExists: true)   // assembly_qc.nf
mito_circular_script = file("${projectDir}/bin/python/plot_mito_circular.py",              checkIfExists: true)   // organelle_assembly.nf
```

For **every** module in the table, the input-block change is the same: add `path(<var>)` as the last `input:` line (e.g. `path(summary_report_script)`, `path(dotplot_script)`, `path(compartments_script)`, …).

---

## Step 4 — Validate

1. Apply, then `nextflow run … -resume`. **Expect these 10 processes to re-run once** — their hash legitimately changed (new declared input + new command string). One-time cost; they're all terminal report/plot steps.
2. **Prove the footgun is gone:** add a trivial comment line to e.g. `bin/r_scripts/compile_qc.R`, then `-resume` again. `COMPILE_FINAL_QC` should now **re-run** (before Track 1 it would have been skipped). That's the whole point.
3. Spot-check a work dir: `.command.sh` shows `Rscript compile_qc.R …` (staged basename, no `${projectDir}` path), and `compile_qc.R` is staged (symlinked) into the task directory.

---

## Notes / gotchas

- **Verify against live files first.** The module input-blocks + invocation lines above are read off the Phase-1 `/mnt/project` snapshot. The invocation lines are stable, but confirm no *new* input was added to a module since then before you add the script input (so the positional order stays right). Call sites are from your current `main.nf`.
- **Leave the inline heredocs alone.** `Rscript - <<'RSCRIPT'` / `python3 <<'PYEOF'` (in `hic_coverage.nf`, `hic_pair_stats.nf`, `teloclip.nf`, `tidk.nf`, `fcs_blob_evidence_report.nf`) and the `Rscript -e '…'` at `pairwise_alignment.nf:41` are embedded in the command string → **already hashed**. Only the 10 external `.R`/`.py` files need staging.
- **`combine_hic_mapping_qc.nf:149`** — the `# Rscript …combine_hic_mapping_qc.R` line is commented out; ignore it.
- The other Track-1 items are separate small jobs (do alongside or after): **dead-include removal** + the **`GAP_FILLING` unquoted-heredoc** fix (§9), and surfacing the **genome-size estimate into `SUMMARY_REPORT`** (the one Phase-4a loose end) — folded into the `SUMMARY_REPORT` touch you're already making here.
