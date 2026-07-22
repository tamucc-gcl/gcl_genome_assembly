#!/usr/bin/env Rscript

#' =============================================================================
#' GENERATE PIPELINE SUMMARY REPORT
#' =============================================================================
#' Reads data passed natively via Nextflow channels (no directory scanning).
#' Produces a Markdown report with:
#'   - Links to final genome assemblies, compiled QC CSV, interactive HTML report
#'   - Assembly QC tables (per-haplotype, not averaged)
#'   - Per-sample visual comparison table (snails, contact maps, within-sample dotplot)
#'   - QC trend plots from COMPILE_FINAL_QC (collapsible)
#'   - Cross-stage comparison table (collapsible)
#'   - Telomere and pairwise alignment summaries (collapsible)
#'
#' The report publishes to {outdir}/ directly, so all relative paths are
#' simply subdir/filename.
#' =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(argparse)
})

# =============================================================================
# CLI arguments
# =============================================================================
parser <- ArgumentParser(description = "Generate pipeline summary report from Nextflow channel data")
parser$add_argument("--manifest",         required = TRUE, help = "Report manifest TSV")
parser$add_argument("--compiled_qc",      required = TRUE, help = "Compiled QC CSV (assembly_qc_metrics.csv)")
parser$add_argument("--telomere_summary", required = TRUE, help = "Telomere summary TSV (or NO_TELOMERES)")
parser$add_argument("--pairwise_summary", required = TRUE, help = "Pairwise alignment summary TSV (or NO_PAIRWISE)")
parser$add_argument("--output",           default = "assembly_report.md", help = "Output markdown file")
parser$add_argument("--img_width",        default = 500, type = "integer", help = "Image display width in pixels")
parser$add_argument("--mito_stats",       required = TRUE, help = "Mitogenome stats TSV (or NO_MITO_STATS)")
parser$add_argument("--teloclip_stats",   required = TRUE, help = "Teloclip extension stats TSV (or NO_TELOCLIP)")
parser$add_argument("--sample_taxonomy", default = "NO_TAXONOMY",    help = "Per-sample taxonomy TSV (sample/taxid/species/kingdom/busco_lineage) or NO_TAXONOMY")
parser$add_argument("--genome_size",     default = "NO_GENOME_SIZE", help = "Per-sample genome-size estimate TSV (sample/est_genome_size_bp) or NO_GENOME_SIZE")
parser$add_argument("--workflow_info", default = "NO_WORKFLOW_INFO", help = "Workflow provenance TSV (key/value) or NO_WORKFLOW_INFO")
parser$add_argument("--run_info",      default = "NO_RUN_INFO",      help = "Per-sample run-info TSV (sample/evidence/strategy) or NO_RUN_INFO")
parser$add_argument("--flag_busco", default = 90, type = "double", help = "Status flag: warn if BUSCO complete percent below this")
parser$add_argument("--flag_qv",    default = 40, type = "double", help = "Status flag: warn if Merqury QV below this")
parser$add_argument("--flag_kmer",  default = 90, type = "double", help = "Status flag: warn if k-mer completeness percent below this")
parser$add_argument("--flag_size_pct", default = 0, type = "double", help = "Status flag: warn if assembled length deviates from estimate by more than this percent (0 = off)")
parser$add_argument("--busco_fallback", default = "eukaryota_odb10", help = "Configured BUSCO fallback lineage (params.busco_lineage), for provenance flagging")
parser$add_argument("--ran_purge_dups", default = "false", help = "Whether purge_dups ran (params.run_purge_dups)")
parser$add_argument("--ran_decontam",   default = "false", help = "Whether FCS decontamination ran (params.decon.run_on_contigs)")
parser$add_argument("--versions", default = "NO_VERSIONS", help = "Software versions TSV (tool/version) or NO_VERSIONS")

args <- parser$parse_args()

# =============================================================================
# Helpers
# =============================================================================

make_markdown_table <- function(df) {
  if (is.null(df) || nrow(df) == 0) return("*No data available.*")
  header <- paste0("| ", paste(colnames(df), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  body <- apply(df, 1, function(row) {
    paste0("| ", paste(row, collapse = " | "), " |")
  })
  c(header, separator, body)
}

make_collapsible <- function(lines, summary_text) {
  c("<details>",
    sprintf("<summary><b>%s</b></summary>", summary_text),
    "", lines, "",
    "</details>", "")
}

rel_path <- function(subdir, filename) {
  file.path(subdir, filename)
}

img_tag <- function(src, alt, width = args$img_width) {
  sprintf('<img src="%s" alt="%s" width="%d">', src, alt, width)
}

# =============================================================================
# Read manifest
# =============================================================================
message("Reading manifest: ", args$manifest)
manifest <- read_tsv(args$manifest, col_types = cols(.default = "c"))

assemblies     <- manifest %>% filter(type == "assembly")
snail_plots    <- manifest %>% filter(type == "snail")
contact_maps   <- manifest %>% filter(type == "contact_map")
dotplots       <- manifest %>% filter(type == "dotplot")
riparian_plots <- manifest %>% filter(type == "riparian")
qc_plots       <- manifest %>% filter(type == "qc_plot")
compiled_csvs  <- manifest %>% filter(type == "compiled_qc")
report_htmls   <- manifest %>% filter(type == "assembly_report_html")
mito_assemblies <- manifest %>% filter(type == "mitogenome")
mito_gene_maps  <- manifest %>% filter(type == "mito_gene_map")
mito_stats_rows <- manifest %>% filter(type == "mito_stats")
genomescope_plots <- manifest %>% filter(type == "genomescope_plot")

all_hap_ids <- assemblies$id %>% sort()
all_sample_ids <- all_hap_ids %>%
  str_replace("_(hap[12]|primary)$", "") %>%
  unique() %>%
  sort()

message(sprintf("Samples: %s", paste(all_sample_ids, collapse = ", ")))

# =============================================================================
# Read compiled QC data
# =============================================================================
message("Reading compiled QC: ", args$compiled_qc)

qc_data <- tryCatch(
  read_csv(args$compiled_qc, show_col_types = FALSE),
  error = function(e) { message("WARNING: ", e$message); NULL }
)

# --- Run Summary (provenance + inputs + stages) ---
wf <- NULL
if (!str_detect(basename(args$workflow_info), "NO_WORKFLOW_INFO") &&
    file.exists(args$workflow_info) && file.size(args$workflow_info) > 0) {
  wf <- tryCatch(
    { d <- read_tsv(args$workflow_info, col_types = cols(.default = "c")); setNames(d$value, d$key) },
    error = function(e) NULL)
}
ri <- NULL
if (!str_detect(basename(args$run_info), "NO_RUN_INFO") &&
    file.exists(args$run_info) && file.size(args$run_info) > 0) {
  ri <- tryCatch(read_tsv(args$run_info, col_types = cols(.default = "c")), error = function(e) NULL)
}

prov <- character()
if (!is.null(wf)) {
  gv <- function(k) if (k %in% names(wf) && !is.na(wf[[k]])) wf[[k]] else "unknown"
  rev_bit <- if (gv("commit") != "unknown")
               sprintf("%s / `%s`", gv("revision"), substr(gv("commit"), 1, 10)) else gv("revision")
  prov <- c(
    sprintf("- **Pipeline:** %s %s (%s)", gv("pipeline"), gv("version"), rev_bit),
    sprintf("- **Nextflow:** %s  ·  **Profile:** %s  ·  **Run:** %s  ·  **Started:** %s",
            gv("nextflow"), gv("profile"), gv("run_name"), gv("start"))
  )
}

n_samples   <- length(all_sample_ids)
sample_line <- sprintf("- **Samples:** %d", n_samples)
if (!is.null(ri) && nrow(ri) > 0) {
  as_lgl <- function(x) tolower(as.character(x)) %in% c("true","t","1","yes")
  parts <- c(sprintf("HiFi %d", sum(as_lgl(ri$hifi))),
             sprintf("Hi-C %d", sum(as_lgl(ri$hic))),
             sprintf("short-read %d", sum(as_lgl(ri$shortread))))
  if (sum(as_lgl(ri$tellseq)) > 0) parts <- c(parts, sprintf("TellSeq %d", sum(as_lgl(ri$tellseq))))
  asm <- ri %>% count(assembler) %>% mutate(l = sprintf("%s (%d)", assembler, n))
  sample_line <- c(
    sprintf("- **Samples:** %d  ·  **Inputs (samples with):** %s", n_samples, paste(parts, collapse = " · ")),
    sprintf("- **Assemblers:** %s", paste(asm$l, collapse = ", "))
  )
}

stages_line <- character()
if (!is.null(qc_data) && nrow(qc_data) > 0 && "stage" %in% names(qc_data)) {
  ord <- c("ctg.base","ctg.mito","ctg.purged","ctg.cor","ctg.deco",
           "scaf.r1","scaf.r2","gap_fill","teloclip","final")
  st  <- unique(qc_data$stage)
  st  <- c(intersect(ord, st), sort(setdiff(st, ord)))
  if (length(st) > 0) stages_line <- sprintf("- **QC stages present:** %s", paste(st, collapse = " → "))
}

run_summary <- c("## Run Summary", "", prov, sample_line, stages_line, "")

# =============================================================================
# Build the report
# =============================================================================
md <- character()

md <- c(md,
        "# Genome Assembly Pipeline Report",
        "",
        sprintf("> Generated: **%s**", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "",
        run_summary,
        "---",
        "",
        "## Table of Contents",
        "",
        "1. [Final Genome Assemblies](#1-final-genome-assemblies)",
        "2. [Assembly QC Summary](#2-assembly-qc-summary)",
        "3. [Visual Summary](#3-visual-summary)",
        "4. [Assembly QC Across Pipeline Stages](#4-assembly-qc-across-pipeline-stages)",
        "5. [Mitochondrial Genome](#5-mitochondrial-genome)",
        "6. [Telomere Detection](#6-telomere-detection)",
        "7. [Pairwise Alignment Summary](#7-pairwise-alignment-summary)",
        "8. [Methods and Citations](#8-methods-and-citations)",
        "",
        "---",
        ""
)

# =============================================================================
# Section 1: Final Genome Assemblies + key resource links
# =============================================================================
md <- c(md,
        "## 1. Final Genome Assemblies",
        "",
        "Final gap-filled, chromosome-level assemblies produced by the [gcl_genome_assembly](https://github.com/tamucc-gcl/gcl_genome_assembly) pipeline.",
        ""
)

# Key resource links
resource_links <- character()
if (nrow(compiled_csvs) > 0) {
  e <- compiled_csvs[1, ]
  resource_links <- c(resource_links,
                      sprintf("- **Compiled QC data (CSV):** [%s](%s)", e$filename, rel_path(e$subdir, e$filename)))
}
if (nrow(report_htmls) > 0) {
  e <- report_htmls[1, ]
  resource_links <- c(resource_links,
                      sprintf("- **Interactive QC report (HTML):** [%s](%s) *(download and open in browser)*",
                              e$filename, rel_path(e$subdir, e$filename)))
}
if (length(resource_links) > 0) {
  md <- c(md, "### Key Resources", "", resource_links, "")
}

# Assembly table
if (nrow(assemblies) > 0) {
  asm_table <- assemblies %>%
    arrange(id) %>%
    mutate(
      sample_id = str_replace(id, "_(hap[12]|primary)$", ""),
      hap = str_extract(id, "hap[12]|primary"),
      hap = if_else(hap == "primary", "hap1", hap),
      link = sprintf("[%s](%s)", filename, rel_path(subdir, filename))
    ) %>%
    select(sample_id, hap, link) %>%
    pivot_wider(names_from = hap, values_from = link, values_fill = "—")
  # Haploid samples only produce a hap1 column; ensure both exist before renaming.
  if (!"hap1" %in% names(asm_table)) asm_table$hap1 <- "—"
  if (!"hap2" %in% names(asm_table)) asm_table$hap2 <- "—"
  asm_table <- asm_table %>%
    select(Sample = sample_id, `Haplotype 1` = hap1, `Haplotype 2` = hap2)
  
  md <- c(md, "### Assembly Files", "", make_markdown_table(asm_table), "")
} else {
  md <- c(md, "*No final assemblies found in manifest.*", "")
}

# ---- Sample taxonomy & genome-size estimate (4b-i) ----
read_side_tsv <- function(path, sentinel) {
  if (str_detect(basename(path), sentinel) || !file.exists(path) || file.size(path) == 0) return(NULL)
  tryCatch(read_tsv(path, col_types = cols(.default = "c")),
           error = function(e) { message("WARNING: ", e$message); NULL })
}

tax_tbl  <- read_side_tsv(args$sample_taxonomy, "NO_TAXONOMY")
size_tbl <- read_side_tsv(args$genome_size,     "NO_GENOME_SIZE")

if (!is.null(tax_tbl) || !is.null(size_tbl)) {
  info <- tibble(sample = all_sample_ids)
  if (!is.null(tax_tbl))  info <- left_join(info, tax_tbl,  by = "sample")
  if (!is.null(size_tbl)) info <- left_join(info, size_tbl, by = "sample")

  if ("est_genome_size_bp" %in% names(info)) {
    info <- info %>%
      mutate(est_genome_size = {
        v <- suppressWarnings(as.numeric(str_remove_all(est_genome_size_bp, "[, ]")))
        ifelse(is.na(v), "—", comma(v, accuracy = 0.01, scale = 1/1e9, suffix = " Gb"))
      }) %>%
      select(-est_genome_size_bp)
  }

  col_map <- c(sample = "Sample", species = "Species", taxid = "Taxid",
               kingdom = "Kingdom", busco_lineage = "BUSCO Lineage",
               est_genome_size = "Est. Genome Size")
  present <- intersect(names(col_map), names(info))
  info_disp <- info %>%
    select(all_of(present)) %>%
    mutate(across(everything(), ~ replace_na(as.character(.x), "—")))
  names(info_disp) <- col_map[present]

  md <- c(md,
          "### Sample Taxonomy & Genome-Size Estimate", "",
          "Per-sample organism identity (resolved from the NCBI taxid) and the k-mer-based",
          "haploid genome-size estimate from [GenomeScope2](https://github.com/tbenavi1/genomescope2.0).", "",
          make_markdown_table(info_disp), "")
}

# =============================================================================
# Section 2: Assembly QC Summary — per-haplotype, not averaged
# =============================================================================
md <- c(md, "## 2. Assembly QC Summary", "")

if (!is.null(qc_data) && nrow(qc_data) > 0) {
  
  # --- Per-haplotype overview table (wide: metrics as rows, samples as columns) ---
  md <- c(md, "### Overview (Final Assemblies)", "")
  
  # First, normalise BUSCO counts to proportions using total_busco
  final_data <- qc_data %>%
    filter(stage == "final") %>%
    mutate(
      across(c(hap1, hap2), ~ case_when(
        metric %in% c("complete", "single", "duplicated", "fragmented", "missing") ~
          . / .[metric == "total_busco"],
        TRUE ~ .
      )),
      .by = sample_id
    )
  
  overview_long <- final_data %>%
    pivot_longer(cols = c(hap1, hap2, both),
                 names_to = "haplotype",
                 values_to = "value",
                 values_drop_na = TRUE) %>%
    filter(metric %in% c("coverage", "Total length", "L90", "auN",
                         "complete", "Largest contig", "GC (%)", "qv",
                         "kmer_completeness"))
  
  if (nrow(overview_long) > 0) {
    # Format values
    fmt_value <- function(value, metric) {
      case_when(
        metric == "complete"                          ~ percent(value, accuracy = 0.1),
        metric == "coverage"                        ~ comma(value, accuracy = 0.1),
        metric == "Total length"                      ~ comma(value, accuracy = 0.01, scale = 1/1e9, suffix = " Gb"),
        metric %in% c("Largest contig", "auN")        ~ comma(value, accuracy = 0.1, scale = 1/1e6, suffix = " Mb"),
        metric %in% c("GC (%)", "kmer_completeness")  ~ percent(value, scale = 1, accuracy = 0.1),
        metric == "qv"                                ~ comma(value, accuracy = 0.1),
        TRUE                                          ~ comma(value, accuracy = 1)
      )
    }
    
    # Metric display names and order
    metric_display <- c(
      "Total length"       = "Total Length",
      "Largest contig"     = "Largest Scaffold",
      "auN"                = "auN",
      "L90"                = "L90",
      "GC (%)"             = "GC (%)",
      "coverage"         = "Coverage",
      "complete"           = "BUSCO Complete",
      "qv"                 = "QV",
      "kmer_completeness"  = "K-mer Completeness"
    )
    metric_order <- names(metric_display)
    
    # Diploid-level metrics (single value per sample)
    diploid_metrics <- c("qv", "kmer_completeness")
    
    # Per-haplotype metrics: format as "hap1 / hap2" slash notation
    # Per-haplotype metrics: format as "hap1 / hap2" slash notation
    per_hap <- overview_long %>%
      filter(!metric %in% diploid_metrics, haplotype %in% c("hap1", "hap2")) %>%
      mutate(value_fmt = fmt_value(value, metric)) %>%
      select(sample_id, metric, haplotype, value_fmt) %>%
      pivot_wider(names_from = haplotype, values_from = value_fmt, values_fill = "—")
    # Haploid samples only produce hap1; add hap2 so the slash notation renders "value / —".
    if (!"hap1" %in% names(per_hap)) per_hap$hap1 <- "—"
    if (!"hap2" %in% names(per_hap)) per_hap$hap2 <- "—"
    per_hap <- per_hap %>%
      mutate(combined = paste(hap1, "/", hap2)) %>%
      select(sample_id, metric, combined)
    
    # Diploid-level metrics: prefer the combined ('both') value; for haploid
    # (no 'both' line) fall back to the single hap1 value.
    diploid <- overview_long %>%
      filter(metric %in% diploid_metrics, haplotype %in% c("both", "hap1", "hap2")) %>%
      mutate(haplotype = factor(haplotype, levels = c("both", "hap1", "hap2"))) %>%
      group_by(sample_id, metric) %>%
      arrange(haplotype, .by_group = TRUE) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(combined = fmt_value(value, metric)) %>%
      select(sample_id, metric, combined)
    
    # Combine and pivot to wide: rows = samples, columns = metrics
    overview_wide <- bind_rows(per_hap, diploid) %>%
      mutate(
        metric = factor(metric, levels = metric_order),
        metric_label = metric_display[as.character(metric)]
      ) %>%
      arrange(sample_id, metric) %>%
      select(sample_id, metric_label, combined) %>%
      pivot_wider(names_from = metric_label, values_from = combined, values_fill = "") %>%
      rename(Sample = sample_id)
    
    # --- Est. Genome Size + % of Estimate (report-improvement batch) ---
    # Est. Genome Size: per-sample haploid estimate (GenomeScope2, via --genome_size).
    # % of Estimate: per-hap (each haplotype's assembled Total length / the sample's
    # haploid estimate) so it reads directly against the "Total Length" row above.
    gs_path <- args$genome_size
    est <- NULL
    if (!str_detect(basename(gs_path), "NO_GENOME_SIZE") &&
        file.exists(gs_path) && file.size(gs_path) > 0) {
      est <- tryCatch(
        read_tsv(gs_path, show_col_types = FALSE) %>%
          transmute(
            sample_id = as.character(sample),
            est_bp    = suppressWarnings(as.numeric(str_remove_all(as.character(est_genome_size_bp), "[, ]")))
          ) %>%
          filter(!is.na(est_bp)),
        error = function(e) { message("WARNING (genome-size cols): ", e$message); NULL }
      )
    }

    if (!is.null(est) && nrow(est) > 0) {
      pct_tbl <- overview_long %>%
        filter(metric == "Total length", haplotype %in% c("hap1", "hap2")) %>%
        select(sample_id, haplotype, total_len = value) %>%
        inner_join(est, by = "sample_id") %>%
        mutate(pct = percent(total_len / est_bp, accuracy = 1)) %>%
        select(sample_id, haplotype, pct) %>%
        pivot_wider(names_from = haplotype, values_from = pct, values_fill = "—")
      if (!"hap1" %in% names(pct_tbl)) pct_tbl$hap1 <- "—"
      if (!"hap2" %in% names(pct_tbl)) pct_tbl$hap2 <- "—"
      pct_tbl <- pct_tbl %>% transmute(sample_id, `% of Estimate` = paste(hap1, "/", hap2))

      size_cols <- est %>%
        transmute(sample_id,
                  `Est. Genome Size` = comma(est_bp, accuracy = 0.01, scale = 1/1e9, suffix = " Gb")) %>%
        left_join(pct_tbl, by = "sample_id")

      overview_wide <- overview_wide %>%
        left_join(size_cols, by = c("Sample" = "sample_id")) %>%
        mutate(across(c(`Est. Genome Size`, `% of Estimate`),
                      ~ replace_na(as.character(.x), "—")))
      if ("Total Length" %in% names(overview_wide)) {
        overview_wide <- overview_wide %>%
          relocate(`Est. Genome Size`, `% of Estimate`, .after = `Total Length`)
      }
    }

    # --- Per-sample QC Status flag (report-improvement batch) ---
    # ✅ = all checks pass · ⚠️ = one or more below the configured threshold. Worst-case
    # (min) across haplotypes per sample. Thresholds are organism-dependent, so they're
    # per-run tunable via --flag_busco / --flag_qv / --flag_kmer / --flag_size_pct.
    thr_busco <- suppressWarnings(as.numeric(args$flag_busco))
    thr_qv    <- suppressWarnings(as.numeric(args$flag_qv))
    thr_kmer  <- suppressWarnings(as.numeric(args$flag_kmer))
    thr_size  <- suppressWarnings(as.numeric(args$flag_size_pct))   # ±% window; <= 0 disables

    verdict_src <- overview_long %>%
      filter(metric %in% c("complete", "qv", "kmer_completeness")) %>%
      group_by(sample_id, metric) %>%
      summarise(value = suppressWarnings(min(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric, values_from = value)
    for (m in c("complete", "qv", "kmer_completeness"))
      if (!m %in% names(verdict_src)) verdict_src[[m]] <- NA_real_

    # Size check: per-hap assembled length / haploid estimate; keep the hap furthest from
    # 1.0. Reuses `est` from the est-size block above; off unless --flag_size_pct > 0
    # (estimates are only meaningful once ploidy is correct).
    size_dev <- tibble(sample_id = character(), size_ratio = numeric())
    if (thr_size > 0 && !is.null(est) && nrow(est) > 0) {
      size_dev <- overview_long %>%
        filter(metric == "Total length", haplotype %in% c("hap1", "hap2")) %>%
        select(sample_id, total_len = value) %>%
        inner_join(est, by = "sample_id") %>%
        mutate(ratio = total_len / est_bp) %>%
        group_by(sample_id) %>%
        slice_max(abs(ratio - 1), n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(sample_id, size_ratio = ratio)
    }
    verdict_src <- verdict_src %>% left_join(size_dev, by = "sample_id")
    if (!"size_ratio" %in% names(verdict_src)) verdict_src$size_ratio <- NA_real_

    # 'complete' is a 0–1 fraction (normalised upstream); qv & kmer_completeness are 0–100.
    verdict <- verdict_src %>%
      rowwise() %>%
      mutate(
        flags = paste(c(
          if (!is.na(complete)          && complete * 100 < thr_busco)   sprintf("BUSCO<%g%%", thr_busco),
          if (!is.na(qv)                && qv < thr_qv)                  sprintf("QV<%g", thr_qv),
          if (!is.na(kmer_completeness) && kmer_completeness < thr_kmer) sprintf("k-mer<%g%%", thr_kmer),
          if (thr_size > 0 && !is.na(size_ratio) &&
              abs(size_ratio - 1) * 100 > thr_size)                      sprintf("size %s", percent(size_ratio, accuracy = 1))
        ), collapse = ", "),
        Status = if (nchar(flags) == 0) "✅" else "⚠️"
      ) %>%
      ungroup() %>%
      select(sample_id, Status, flags)

    overview_wide <- overview_wide %>%
      left_join(verdict %>% select(sample_id, Status), by = c("Sample" = "sample_id")) %>%
      mutate(Status = replace_na(Status, "—")) %>%
      relocate(Status, .after = Sample)

    md <- c(md, make_markdown_table(overview_wide), "")

    crit <- sprintf("BUSCO complete ≥ %g%%, QV ≥ %g, k-mer completeness ≥ %g%%",
                    thr_busco, thr_qv, thr_kmer)
    if (thr_size > 0)
      crit <- paste0(crit, sprintf(", assembled length within ±%g%% of estimate", thr_size))
    md <- c(md,
      paste0("> **Status:** ✅ all checks pass · ⚠️ below threshold (", crit,
             "). Thresholds are per-run tunable."),
      "")
    flagged <- verdict %>% filter(nchar(flags) > 0)
    if (nrow(flagged) > 0) {
      md <- c(md,
        paste0("Flagged: ",
               paste(sprintf("**%s** (%s)", flagged$sample_id, flagged$flags), collapse = "; ")),
        "")
    }

    flagged <- verdict %>% filter(nchar(flags) > 0)
    if (nrow(flagged) > 0) {
      md <- c(md,
        paste0("Flagged: ",
               paste(sprintf("**%s** (%s)", flagged$sample_id, flagged$flags), collapse = "; ")),
        "")
    }

    # --- BUSCO lineage provenance (report-improvement batch) ---
    # Which lineage each BUSCO % was scored against, and how it was chosen: per-sample
    # from the NCBI taxid (section 1), else the configured fallback.
    if (!is.null(tax_tbl) && "busco_lineage" %in% names(tax_tbl)) {
      lin <- tax_tbl %>%
        filter(sample %in% overview_wide$Sample) %>%
        transmute(sample,
                  lineage = ifelse(is.na(busco_lineage) | busco_lineage == "", "—", busco_lineage))
      fb <- tolower(trimws(as.character(args$busco_fallback)))
      uniq_lin <- sort(unique(lin$lineage[lin$lineage != "—"]))

      note <- paste0("> **BUSCO lineage:** completeness is scored per sample against the lineage ",
                     "resolved from its NCBI taxid (section 1); where a specific clade can't be ",
                     "resolved it falls back to `params.busco_lineage`")
      if (nchar(fb) > 0) note <- paste0(note, sprintf(" (here `%s`)", args$busco_fallback))
      note <- paste0(note, ".")
      if (length(uniq_lin) == 1) {
        note <- paste0(note, sprintf(" All samples scored against `%s`.", uniq_lin))
      } else if (length(uniq_lin) > 1) {
        pairs <- lin %>% arrange(sample) %>%
          mutate(s = sprintf("%s → `%s`", sample, lineage)) %>% pull(s)
        note <- paste0(note, " Per sample: ", paste(pairs, collapse = "; "), ".")
      }
      md <- c(md, note, "")

      if (nchar(fb) > 0) {
        on_fb <- sort(lin$sample[tolower(lin$lineage) == fb])
        if (length(on_fb) > 0) {
          md <- c(md,
            sprintf("⚠️ Scored against the broad fallback lineage (completeness may read high): %s.",
                    paste(sprintf("**%s**", on_fb), collapse = ", ")),
            "")
        }
      }
    }
  }
  
  # --- Detailed final assembly table (collapsible) ---
  final_detail <- final_data %>%
    pivot_longer(cols = c(hap1, hap2, both),
                 names_to = "haplotype",
                 values_to = "value",
                 values_drop_na = TRUE) %>%
    filter(!(metric %in% c("error_rate", "error_kmer", "total_assembly_kmer") & analysis == "merqury"),
           !(metric %in% c("kmer_found", "total_hifi_kmer") & analysis == "merqury"),
           metric != "total_busco",
           !str_detect(metric, ">=")) %>%
    mutate(value_fmt = case_when(
      analysis == "busco" ~ percent(value, accuracy = 0.1),
      metric == "coverage" ~ comma(value, accuracy = 0.1),
      metric %in% c("GC (%)", "kmer_completeness") ~ percent(value, accuracy = 0.1, scale = 1),
      TRUE ~ comma(value, accuracy = 1)
    )) %>%
    select(sample_id, metric, haplotype, value_fmt) %>%
    pivot_wider(names_from = haplotype, values_from = value_fmt, values_fill = "") %>%
    mutate(sample_id = if_else(lag(sample_id, default = "") == sample_id, "", sample_id)) %>%
    rename(Sample = sample_id, Metric = metric)
  
  if (nrow(final_detail) > 0) {
    md <- c(md,
            make_collapsible(
              make_markdown_table(final_detail),
              "Click to expand: All final assembly metrics"
            )
    )
  }
}

# =============================================================================
# Section 3: Visual Summary — per-sample table with hap1/hap2 columns
# =============================================================================
md <- c(md,
        "## 3. Visual Summary",
        "",
        "Each row shows one sample with snail plots, contact maps, the within-sample",
        "haplotype-vs-haplotype dotplot, and riparian (ribbon) synteny plot.",
        ""
)

# We build an HTML table because markdown tables can't embed images well at
# controlled sizes, and we need multi-column image layouts.

# Pick the best contact map resolution to display (largest resolution = most zoomed out)
# The filename pattern is {hap_id}_{stage}_{resolution}bp_contact_map.png
# Pick the one with the largest number (e.g., 1000000bp)
pick_best_contact_map <- function(cmaps_for_hap) {
  if (nrow(cmaps_for_hap) == 0) return(NULL)
  cmaps_for_hap %>%
    mutate(res = str_extract(filename, "[0-9]+(?=bp_contact_map)") %>% as.numeric()) %>%
    arrange(desc(res)) %>%
    slice(1)
}

has_snails <- nrow(snail_plots) > 0
has_cmaps  <- nrow(contact_maps) > 0
has_dots   <- nrow(dotplots) > 0
has_riparian <- nrow(riparian_plots) > 0

# Determine which columns to include
col_headers <- c("Sample")
if (has_snails) col_headers <- c(col_headers, "Hap1 Snail", "Hap2 Snail")
if (has_cmaps)  col_headers <- c(col_headers, "Hap1 Contact Map", "Hap2 Contact Map")
if (has_dots)   col_headers <- c(col_headers, "Hap1 vs Hap2 Dotplot")
if (has_riparian) col_headers <- c(col_headers, "Hap1 vs Hap2 Riparian")

if (has_snails || has_cmaps || has_dots || has_riparian) {
  
  img_w <- args$img_width
  
  # Start HTML table
  md <- c(md, "<table>", "<tr>")
  for (h in col_headers) {
    md <- c(md, sprintf("  <th>%s</th>", h))
  }
  md <- c(md, "</tr>")
  
  for (sid in all_sample_ids) {
    # Resolve this sample's haplotype ids from the manifest instead of assuming
    # hap1/hap2: diploid -> [sid_hap1, sid_hap2], haploid -> [sid_primary].
    sample_haps <- sort(all_hap_ids[str_replace(all_hap_ids, "_(hap[12]|primary)$", "") == sid])
    slot1_id <- sample_haps[1]
    slot2_id <- if (length(sample_haps) >= 2) sample_haps[2] else NA_character_
    
    md <- c(md, "<tr>")
    
    # Sample name
    md <- c(md, sprintf('  <td><b>%s</b></td>', sid))
    
    # Snail plots (slot 1 = hap1/primary, slot 2 = hap2 or — for haploid)
    if (has_snails) {
      for (hid in c(slot1_id, slot2_id)) {
        snail_row <- if (!is.na(hid)) snail_plots %>% filter(id == hid) else snail_plots[0, ]
        if (nrow(snail_row) > 0) {
          src <- rel_path(snail_row$subdir[1], snail_row$filename[1])
          md <- c(md, sprintf('  <td>%s</td>', img_tag(src, hid, width = img_w)))
        } else {
          md <- c(md, "  <td>—</td>")
        }
      }
    }
    
    # Contact maps (pick best resolution)
    if (has_cmaps) {
      for (hid in c(slot1_id, slot2_id)) {
        hap_cmaps <- if (!is.na(hid)) contact_maps %>% filter(id == hid) else contact_maps[0, ]
        best <- pick_best_contact_map(hap_cmaps)
        if (!is.null(best) && nrow(best) > 0) {
          src <- rel_path(best$subdir[1], best$filename[1])
          md <- c(md, sprintf('  <td>%s</td>', img_tag(src, hid, width = img_w)))
        } else {
          md <- c(md, "  <td>—</td>")
        }
      }
    }
    
    # Within-sample dotplot (hap1 vs hap2) — diploid only; haploid has no pair
    if (has_dots) {
      dot_row <- if (!is.na(slot2_id)) {
        dotplots %>%
          filter(
            (id == slot1_id & id2 == slot2_id) |
              (id == slot2_id & id2 == slot1_id)
          )
      } else dotplots[0, ]
      if (nrow(dot_row) > 0) {
        src <- rel_path(dot_row$subdir[1], dot_row$filename[1])
        md <- c(md, sprintf('  <td>%s</td>', img_tag(src, paste(slot1_id, "vs", slot2_id), width = img_w)))
      } else {
        md <- c(md, "  <td>—</td>")
      }
    }
    
    # Within-sample riparian plot (hap1 vs hap2) — diploid only
    if (has_riparian) {
      rip_row <- if (!is.na(slot2_id)) {
        riparian_plots %>%
          filter(
            (id == slot1_id & id2 == slot2_id) |
              (id == slot2_id & id2 == slot1_id)
          )
      } else riparian_plots[0, ]
      if (nrow(rip_row) > 0) {
        src <- rel_path(rip_row$subdir[1], rip_row$filename[1])
        md <- c(md, sprintf('  <td>%s</td>', img_tag(src, paste(slot1_id, "vs", slot2_id, "riparian"), width = img_w)))
      } else {
        md <- c(md, "  <td>—</td>")
      }
    }
    
    md <- c(md, "</tr>")
  }
  
  md <- c(md, "</table>", "")
  
  # Link to contact maps folder for all resolutions
  if (has_cmaps) {
    cmap_subdir <- contact_maps$subdir[1]
    md <- c(md,
            sprintf("All contact map resolutions: [%s/](%s/)", cmap_subdir, cmap_subdir),
            "")
  }
  
  # Link to pairwise alignments folder for cross-sample dotplots
  if (has_dots) {
    dot_subdir <- dotplots$subdir[1]
    n_total <- nrow(dotplots)
    n_within <- sum(
      dotplots %>%
        mutate(s1 = str_replace(id, "_(hap[12]|primary)$", ""),
               s2 = str_replace(id2, "_(hap[12]|primary)$", "")) %>%
        pull(s1) == dotplots %>%
        mutate(s2 = str_replace(id2, "_(hap[12]|primary)$", "")) %>%
        pull(s2)
    )
    n_cross <- n_total - n_within
    
    if (n_cross > 0) {
      md <- c(md,
              sprintf("Additional cross-sample dotplots (%d): [%s/](%s/)",
                      n_cross, dot_subdir, dot_subdir),
              "")
    }
  }
}

# --- GenomeScope profiles (per-sample k-mer spectra + model fit) ---
if (nrow(genomescope_plots) > 0) {
  gs_w <- args$img_width
  md <- c(md,
    "### GenomeScope Profiles", "",
    "K-mer spectra with the fitted [GenomeScope2](https://github.com/tbenavi1/genomescope2.0)",
    "model — visual support for the haploid genome-size estimate (peak position), plus",
    "heterozygosity and repeat content. One profile per sample.", "")

  gs <- genomescope_plots %>% arrange(id)
  md <- c(md, "<table>", "<tr><th>Sample</th><th>GenomeScope Profile</th></tr>")
  for (i in seq_len(nrow(gs))) {
    src <- rel_path(gs$subdir[i], gs$filename[i])
    md <- c(md, sprintf("<tr><td><b>%s</b></td><td>%s</td></tr>",
                        gs$id[i], img_tag(src, gs$id[i], width = gs_w)))
  }
  md <- c(md, "</table>", "")
}

# =============================================================================
# Section 4: QC trend plots + cross-stage table (only meaningful with >1 stage)
# =============================================================================
n_stages <- if (!is.null(qc_data) && nrow(qc_data) > 0 && "stage" %in% names(qc_data))
              n_distinct(qc_data$stage) else 0L

if (!is.null(qc_data) && nrow(qc_data) > 0) {
  md <- c(md, "## 4. Assembly QC Across Pipeline Stages", "")

  if (n_stages > 1) {
    plot_lines <- character()

    # Trend plots only if COMPILE_FINAL_QC emitted them
    if (nrow(qc_plots) > 0) {
      plot_lines <- c(plot_lines, "### QC Trend Plots", "")
      plot_order <- c("busco", "kmer", "quast_misc", "contig_count", "contig_length", "trans_cis")
      plot_labels <- c(
        busco = "BUSCO Completeness",
        kmer = "K-mer QV & Completeness",
        quast_misc = "QUAST Assembly Metrics",
        contig_count = "Contig/Scaffold Counts",
        contig_length = "Assembly Size",
        trans_cis = "Hi-C Trans:Cis Ratio"
      )
      ordered_plots <- qc_plots %>%
        mutate(plot_key = str_remove(filename, "\\.png$")) %>%
        mutate(order = match(plot_key, plot_order)) %>%
        arrange(order, filename)
      for (i in seq_len(nrow(ordered_plots))) {
        row <- ordered_plots[i, ]
        src <- rel_path(row$subdir, row$filename)
        label <- plot_labels[row$plot_key]
        if (is.na(label)) label <- row$plot_key
        plot_lines <- c(plot_lines,
                        sprintf("#### %s", label), "",
                        img_tag(src, label, width = 800), "")
      }
    }

    # Cross-stage comparison table (normalise BUSCO before formatting)
    cross_stage <- qc_data %>%
      mutate(
        across(c(hap1, hap2), ~ case_when(
          metric %in% c("complete", "single", "duplicated", "fragmented", "missing") ~
            . / .[metric == "total_busco"],
          TRUE ~ .
        )),
        .by = c(sample_id, stage)
      ) %>%
      pivot_longer(cols = c(hap1, hap2, both),
                   names_to = "haplotype",
                   values_to = "value",
                   values_drop_na = TRUE) %>%
      filter(!(metric %in% c("error_rate", "error_kmer", "total_assembly_kmer") & analysis == "merqury"),
             !(metric %in% c("kmer_found", "total_hifi_kmer") & analysis == "merqury"),
             metric != "total_busco",
             !str_detect(metric, ">=")) %>%
      mutate(value_fmt = case_when(
        analysis == "busco" ~ percent(value, accuracy = 0.1),
        metric == "coverage" ~ comma(value, accuracy = 0.1),
        metric %in% c("GC (%)", "kmer_completeness") ~ percent(value, accuracy = 0.1, scale = 1),
        TRUE ~ comma(value, accuracy = 1)
      ))

    if (nrow(cross_stage) > 0) {
      cross_table <- cross_stage %>%
        select(sample_id, metric, stage, haplotype, value_fmt) %>%
        pivot_wider(names_from = c(stage, haplotype),
                    values_from = value_fmt, values_fill = "") %>%
        select(sample_id, metric,
               matches("_hap1$"), matches("_hap2$"), matches("_both$")) %>%
        mutate(sample_id = if_else(lag(sample_id, default = "") == sample_id, "", sample_id)) %>%
        rename(Sample = sample_id, Metric = metric)

      plot_lines <- c(plot_lines,
                      make_collapsible(
                        c("#### Cross-Stage Metrics Table", "",
                          make_markdown_table(cross_table)),
                        "Click to expand: Assembly metrics across all pipeline stages"
                      )
      )
    }

    md <- c(md,
            make_collapsible(plot_lines,
                             "Click to expand: QC trend plots and cross-stage comparison")
    )
  } else {
    md <- c(md,
      paste("*Only the final QC stage is present, so per-stage trend plots and",
            "cross-stage comparisons aren't shown — they populate automatically once",
            "earlier stages are recorded (purge_dups, scaffolding, gap-filling, etc.).",
            "Full final-assembly metrics are in [section 2](#2-assembly-qc-summary).*"),
      "")
  }
}

# =============================================================================
# Section 5: Mitochondrial Genome
# =============================================================================
# This combines the stats table, FASTA links, and gene map images
# into a single HTML table with one row per sample.
# =============================================================================

md <- c(md, "## 5. Mitochondrial Genome", "")

mito_stats <- NULL
if (!grepl("NO_MITO_STATS", args$mito_stats, fixed = TRUE) &&
    file.exists(args$mito_stats)) {
  mito_stats <- tryCatch(
    read_tsv(args$mito_stats, show_col_types = FALSE),
    error = function(e) { message("WARNING: Could not read mito stats: ", e$message); NULL }
  )
}

if (is.null(mito_stats) || nrow(mito_stats) == 0) {
  md <- c(md, "*No mitochondrial genome assembly data available.*", "")
} else {

  # Build a single HTML table: Sample | Gene Map | Length | Circular | Genes | tRNAs | rRNAs | FASTA
  md <- c(md,
    "<table>",
    "<tr>",
    "  <th>Sample</th>",
    "  <th>Gene Map</th>",
    "  <th>Length (bp)</th>",
    "  <th>Circular</th>",
    "  <th>Genes</th>",
    "  <th>tRNAs</th>",
    "  <th>rRNAs</th>",
    "  <th>GenBank</th>",
    "</tr>"
  )

  for (i in seq_len(nrow(mito_stats))) {
    row <- mito_stats[i, ]
    sid <- row$sample_id

    # Gene map image from manifest
    map_row <- mito_gene_maps %>% filter(id == sid)
    map_cell <- if (nrow(map_row) > 0) {
      src <- rel_path(map_row$subdir[1], map_row$filename[1])
      img_tag(src, paste("Mito gene map:", sid), width = args$img_width)
    } else { "—" }

    # FASTA link from manifest
    mito_gb <- manifest %>% filter(type == "mito_genbank", id == sid)
    gb_cell <- if (nrow(mito_gb) > 0) {
      sprintf('<a href="%s">%s</a>',
              rel_path(mito_gb$subdir[1], mito_gb$filename[1]),
              mito_gb$filename[1])
    } else { "—" }

    # Format stats
    length_fmt <- scales::comma(as.numeric(row$mitogenome_length))
    circular_fmt <- ifelse(tolower(row$circular) == "yes", "Yes", row$circular)

    md <- c(md, sprintf(
      "<tr><td><b>%s</b></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
      sid, map_cell, length_fmt, circular_fmt,
      row$gene_count, row$trna_count, row$rrna_count, gb_cell
    ))
  }

  md <- c(md, "</table>", "")

  # Mito filtering note
  md <- c(md,
    "Mitochondrial contigs were identified and removed from each haplotype",
    "assembly prior to purge\\_dups, Inspector, decontamination, and scaffolding.",
    ""
  )
}

# =============================================================================
# Section 6: Telomere Detection
# =============================================================================
# This section includes:
#   a) Teloclip extension summary (if teloclip was run)
#   b) tidk telomere presence/absence summary table
#   c) tidk density SVG plots per haplotype (from manifest)
# =============================================================================

md <- c(md, "## 6. Telomere Detection", "")

# ---- 6a: Teloclip extension summary ----
teloclip_path <- args$teloclip_stats
has_teloclip <- !str_detect(basename(teloclip_path), "NO_TELOCLIP") &&
  file.exists(teloclip_path) && file.size(teloclip_path) > 0

if (has_teloclip) {
  teloclip_data <- tryCatch(read_tsv(teloclip_path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(teloclip_data) && nrow(teloclip_data) > 0) {

    # Assign each extended contig to its haplotype by PREFIX-matching the known
    # hap ids. teloclip prefixes every contig with "{meta.id}_", so this is robust
    # to the contig suffix (hifiasm _ptgNl, YaHS _scaffold_N, SPAdes _NODE_..).
    # (The old code stripped "_scaffold_N", which unscaffolded HiFi contigs don't
    #  have -> every contig became its own group -> one row per contig.)
    assign_hap <- function(ctg) {
      hits <- all_hap_ids[startsWith(ctg, paste0(all_hap_ids, "_"))]
      if (length(hits) == 0) return(NA_character_)
      hits[which.max(nchar(hits))]           # longest match wins
    }

    tc <- teloclip_data %>%
      filter(extension_length > 0) %>%
      mutate(haplotype_id = vapply(contig, assign_hap, character(1))) %>%
      filter(!is.na(haplotype_id))

    teloclip_summary <- tc %>%
      group_by(haplotype_id) %>%
      summarise(
        contigs_extended  = n_distinct(contig),
        extensions        = n(),
        total_bp_added    = sum(extension_length),
        mean_extension_bp = round(mean(extension_length)),
        max_extension_bp  = max(extension_length),
        .groups = "drop"
      ) %>%
      arrange(haplotype_id)

    md <- c(md,
            "#### Telomere Extension (teloclip)", "",
            "Soft-clipped HiFi read overhangs containing telomeric motifs were used to",
            "extend contig/scaffold ends missing telomeric sequence. Per-haplotype summary:", ""
    )

    if (nrow(teloclip_summary) > 0) {
      tc_table <- teloclip_summary %>%
        mutate(
          total_bp_added    = comma(total_bp_added),
          mean_extension_bp = comma(mean_extension_bp),
          max_extension_bp  = comma(max_extension_bp)
        ) %>%
        rename(
          Haplotype             = haplotype_id,
          `Contigs Extended`    = contigs_extended,
          `Total Extensions`    = extensions,
          `Total bp Added`      = total_bp_added,
          `Mean Extension (bp)` = mean_extension_bp,
          `Max Extension (bp)`  = max_extension_bp
        )
      md <- c(md, make_markdown_table(tc_table), "")

      tc_detail <- tc %>%
        transmute(
          Haplotype = haplotype_id,
          Contig    = str_remove(contig, paste0("^", haplotype_id, "_")),
          End       = end,
          `Extension (bp)` = comma(extension_length)
        ) %>%
        arrange(Haplotype, Contig)
      md <- c(md,
              make_collapsible(
                make_markdown_table(tc_detail),
                sprintf("Click to expand: per-contig extension detail (%s extensions)", nrow(tc_detail))
              )
      )
    }
  }
} else {
  md <- c(md,
          "*Teloclip extension was not run. Enable with `--run_teloclip_extend true`.*", ""
  )
}

# ---- 6b: tidk telomere summary table ----
telo_path <- args$telomere_summary
has_telo <- !str_detect(basename(telo_path), "NO_TELOMERES") &&
  file.exists(telo_path) && file.size(telo_path) > 0

if (has_telo) {
  telo_data <- tryCatch(read_tsv(telo_path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(telo_data) && nrow(telo_data) > 0) {
    md <- c(md,
            "#### Telomere Presence Summary (tidk)", "",
            "Telomeric repeat detection across scaffold ends using",
            "[tidk](https://github.com/tolkit/telomeric-identifier).", ""
    )
    
    telo_display <- telo_data %>%
      mutate(
        total_length = scales::comma(as.numeric(total_length)),
        pct_with_telomere = paste0(pct_with_telomere, "%")
      ) %>%
      rename(
        Haplotype = haplotype_id,
        Scaffolds = scaffolds,
        `Total Length` = total_length,
        `5' Telo` = telomere_5prime,
        `3' Telo` = telomere_3prime,
        Both = telomere_both,
        `5' Only` = telomere_5prime_only,
        `3' Only` = telomere_3prime_only,
        None = telomere_none,
        `% with Telomere` = pct_with_telomere
      )
    
    md <- c(md, make_markdown_table(telo_display), "")
    
    # Per-scaffold detail in collapsible
    telo_detail_path <- file.path(dirname(telo_path), "all_telomeres.tsv")
    if (file.exists(telo_detail_path) && file.size(telo_detail_path) > 0) {
      telo_detail <- tryCatch(read_tsv(telo_detail_path, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(telo_detail) && nrow(telo_detail) > 0) {
        # Only show scaffolds that HAVE telomeres (to keep it compact)
        telo_has <- telo_detail %>%
          filter(telomere_5prime == 1 | telomere_3prime == 1) %>%
          mutate(
            length = scales::comma(as.numeric(length)),
            telomere_5prime = ifelse(telomere_5prime == 1, "Yes", "—"),
            telomere_3prime = ifelse(telomere_3prime == 1, "Yes", "—"),
            telomere_both   = ifelse(telomere_both == 1, "Yes", "—")
          ) %>%
          rename(
            Haplotype = haplotype_id, Scaffold = scaffold, Length = length,
            `5'` = telomere_5prime, `3'` = telomere_3prime, Both = telomere_both
          )
        
        if (nrow(telo_has) > 0) {
          md <- c(md,
                  make_collapsible(
                    make_markdown_table(telo_has),
                    "Click to expand: Scaffolds with detected telomeres"
                  )
          )
        }
      }
    }
  }
} else {
  md <- c(md, "*No telomere detection data available.*", "")
}

# ---- 6c: tidk density plots ----
#tidk_plots <- manifest %>% filter(type == "tidk_plot")
#if (nrow(tidk_plots) > 0) {
#  md <- c(md, "#### Telomere Density Plots (tidk)", "")
#  
#  # Build an HTML table with one column per haplotype within each sample
#  tidk_by_sample <- tidk_plots %>%
#    mutate(
#      sample_id = str_replace(id, "_hap[12]$", ""),
#      hap = str_extract(id, "hap[12]$")
#    ) %>%
#    arrange(sample_id, hap)
#  
#  for (sid in unique(tidk_by_sample$sample_id)) {
#    sample_plots <- tidk_by_sample %>% filter(sample_id == sid)
#    md <- c(md, sprintf("**%s**", sid), "")
#    
#    for (i in seq_len(nrow(sample_plots))) {
#      row <- sample_plots[i, ]
#      src <- rel_path(row$subdir, row$filename)
#      alt <- sprintf("tidk telomere density: %s", row$id)
#      md <- c(md,
#              sprintf("*%s*", row$hap),
#              img_tag(src, alt, width = args$img_width),
#              ""
#      )
#    }
#  }
#}

# ---- End of Section 6 ----

# =============================================================================
# Section 7: Pairwise Alignment Summary (collapsible)
# =============================================================================
pw_path <- args$pairwise_summary
if (!str_detect(basename(pw_path), "NO_PAIRWISE") &&
    file.exists(pw_path) && file.size(pw_path) > 0) {
  pw_data <- tryCatch(read_tsv(pw_path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(pw_data) && nrow(pw_data) > 0) {
    md <- c(md,
            "## 7. Pairwise Alignment Summary", "",
            make_collapsible(
              make_markdown_table(pw_data %>% mutate(across(everything(), as.character))),
              "Click to expand: Pairwise alignment metrics"
            )
    )
  }
}

# =============================================================================
# Section 8: Methods and Citations
# =============================================================================
as_lgl <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes")

sig_sr   <- !is.null(ri) && any(as_lgl(ri$shortread))
sig_hifi <- !is.null(ri) && any(as_lgl(ri$hifi))
sig_hic  <- (!is.null(ri) && any(as_lgl(ri$hic))) || nrow(contact_maps) > 0
sig_mito <- nrow(mito_stats_rows) > 0
sig_syn  <- nrow(dotplots) > 0 || nrow(riparian_plots) > 0
has_teloclip <- !str_detect(basename(args$teloclip_stats), "NO_TELOCLIP") &&
                file.exists(args$teloclip_stats) && file.size(args$teloclip_stats) > 0
asm_set   <- if (!is.null(ri)) sort(unique(tolower(ri$assembler))) else character()
ran_purge <- as_lgl(args$ran_purge_dups)
ran_decon <- as_lgl(args$ran_decontam)

asm_pretty <- c(hifiasm = "hifiasm", spades = "SPAdes")
asm_names  <- unname(ifelse(asm_set %in% names(asm_pretty), asm_pretty[asm_set], asm_set))
asm_lbl <- if (length(asm_names) == 0) "the configured assembler" else
           if (length(asm_names) == 1) asm_names else
           paste(paste(asm_names[-length(asm_names)], collapse = ", "), "and", asm_names[length(asm_names)])

# --- Methods narrative (reflects this run) ---
narr <- sprintf("This run processed %d sample%s", n_samples, if (n_samples == 1) "" else "s")
if (!is.null(ri)) {
  inp <- c()
  if (sig_hifi) inp <- c(inp, sprintf("HiFi (%d)", sum(as_lgl(ri$hifi))))
  if (sig_hic)  inp <- c(inp, sprintf("Hi-C (%d)", sum(as_lgl(ri$hic))))
  if (sig_sr)   inp <- c(inp, sprintf("short-read (%d)", sum(as_lgl(ri$shortread))))
  if (!is.null(ri$tellseq) && any(as_lgl(ri$tellseq)))
    inp <- c(inp, sprintf("TellSeq (%d)", sum(as_lgl(ri$tellseq))))
  if (length(inp)) narr <- paste0(narr, sprintf(" (%s)", paste(inp, collapse = ", ")))
}
narr <- paste0(narr, sprintf(". Genome assembly used %s.", asm_lbl))
if (sig_sr)    narr <- paste0(narr, " Short-read contigs were reduced, scaffolded and gap-closed with Redundans.")
if (ran_purge) narr <- paste0(narr, " Haplotypic duplication was removed with purge_dups.")
narr <- paste0(narr, " Haploid genome size and heterozygosity were estimated from k-mer spectra (Jellyfish + GenomeScope 2).")
if (ran_decon) narr <- paste0(narr, " Assemblies were screened for contaminants with NCBI FCS (FCS-Adaptor and FCS-GX).")
if (sig_hic)   narr <- paste0(narr, " Contigs were scaffolded against Hi-C data with YaHS.")
if (sig_mito)  narr <- paste0(narr, " Organelle genomes were assembled with MitoHiFi.")
narr <- paste0(narr, " Assembly quality was assessed with BUSCO (per-sample lineage; see section 2), Merqury (consensus QV and k-mer completeness) and QUAST (contiguity), with read coverage from minimap2/SAMtools alignments; telomeric repeats were surveyed with tidk.")
if (has_teloclip) narr <- paste0(narr, " Scaffold ends were extended into telomeric repeats with teloclip.")
if (sig_syn)   narr <- paste0(narr, " Synteny was visualised from minimap2 alignments (gggenomes).")
narr <- paste0(narr, " Per-step parameters and exact software versions are recorded in the pipeline's Nextflow execution reports.")

# --- Tool reference registry ---
refs <- c(
  pipeline    = "**gcl_genome_assembly** — TAMU-CC Genomics Core Lab. https://github.com/tamucc-gcl/gcl_genome_assembly",
  hifiasm     = "**hifiasm** — Cheng et al. (2021) *Nat. Methods* 18:170-175.",
  spades      = "**SPAdes** — Prjibelski et al. (2020) *Curr. Protoc. Bioinformatics* 70:e102.",
  redundans   = "**Redundans** — Pryszcz & Gabaldon (2016) *Nucleic Acids Res.* 44:e113.",
  purge       = "**purge_dups** — Guan et al. (2020) *Bioinformatics* 36:2896-2898.",
  yahs        = "**YaHS** — Zhou et al. (2023) *Bioinformatics* 39:btac808.",
  fcs         = "**NCBI FCS** (FCS-GX, FCS-Adaptor) — Astashyn et al. (2024) *Genome Biol.* 25:60.",
  mitohifi    = "**MitoHiFi** — Uliano-Silva et al. (2023) *BMC Bioinformatics* 24:288.",
  jellyfish   = "**Jellyfish** — Marçais & Kingsford (2011) *Bioinformatics* 27:764-770.",
  genomescope = "**GenomeScope 2** — Ranallo-Benavidez et al. (2020) *Nat. Commun.* 11:1432.",
  busco       = "**BUSCO** — Simão et al. (2021) *Bioinformatics* 31:3210-3212.",
  merqury     = "**Merqury** — Rhie et al. (2020) *Genome Biol.* 21:245.",
  quast       = "**QUAST** — Gurevich et al. (2013) *Bioinformatics* 29:1072-1075.",
  minimap2    = "**minimap2** — Li (2018) *Bioinformatics* 34:3094-3100.",
  samtools    = "**SAMtools** — Danecek et al. (2021) *GigaScience* 10:giab008.",
  tidk        = "**tidk** — Brown et al. (2021) *Bioinformatics* 41:btaf049.",
  teloclip    = "**teloclip** — Taranto. https://github.com/Adamtaranto/teloclip",
  gggenomes   = "**gggenomes** — Hackl et al. (2024) *arXiv* arXiv:2411.13556"
)

# --- Select citations matching this run ---
keys <- "pipeline"
if (sig_hifi || "hifiasm" %in% asm_set) keys <- c(keys, "hifiasm")
if ("spades" %in% asm_set)              keys <- c(keys, "spades")
if (sig_sr)                             keys <- c(keys, "redundans")
if (ran_purge)                          keys <- c(keys, "purge")
if (ran_decon)                          keys <- c(keys, "fcs")
if (sig_hic)                            keys <- c(keys, "yahs")
if (sig_mito)                           keys <- c(keys, "mitohifi")
keys <- c(keys, "jellyfish", "genomescope", "busco", "merqury", "quast", "minimap2", "samtools", "tidk")
if (has_teloclip)                       keys <- c(keys, "teloclip")
if (sig_syn)                            keys <- c(keys, "gggenomes")
keys <- unique(keys)

md <- c(md,
        "## 8. Methods and Citations", "",
        narr, "",
        "### Tool References", "",
        "Primary references for the tools used in this run. Exact versions are recorded in the pipeline's Nextflow execution reports (not reproduced here).", "",
        paste0("- ", refs[keys]), "")

# --- Software Versions (auto-populates once per-process version capture is enabled) ---
ver_path <- args$versions
vers <- NULL
if (!str_detect(basename(ver_path), "NO_VERSIONS") &&
    file.exists(ver_path) && file.size(ver_path) > 0) {
  vers <- tryCatch(read_tsv(ver_path, col_types = cols(.default = "c")),
                   error = function(e) { message("WARNING (versions): ", e$message); NULL })
}
md <- c(md, "### Software Versions", "")
if (!is.null(vers) && nrow(vers) >= 1 && ncol(vers) >= 2) {
  names(vers)[1:2] <- c("Tool", "Version")
  vers_disp <- vers %>% distinct(Tool, .keep_all = TRUE) %>%
    arrange(tolower(Tool)) %>% select(Tool, Version)
  md <- c(md, make_markdown_table(vers_disp), "")
} else {
  md <- c(md,
    "*Per-tool versions were not recorded for this run; see the pipeline's Nextflow execution report for the resolved software environment.*", "")
}

# =============================================================================
# Footer
# =============================================================================
md <- c(md, "---", "",
        sprintf("*Report generated on %s by the Genome Assembly Pipeline.*",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

writeLines(md, args$output)
message(sprintf("Report written to: %s (%d lines)", args$output, length(md)))