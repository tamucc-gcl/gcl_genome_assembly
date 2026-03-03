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
qc_plots       <- manifest %>% filter(type == "qc_plot")
compiled_csvs  <- manifest %>% filter(type == "compiled_qc")
report_htmls   <- manifest %>% filter(type == "assembly_report_html")

all_hap_ids <- assemblies$id %>% sort()
all_sample_ids <- all_hap_ids %>%
  str_replace("_hap[12]$", "") %>%
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

# =============================================================================
# Build the report
# =============================================================================
md <- character()

md <- c(md,
  "# Genome Assembly Pipeline Report",
  "",
  sprintf("> Generated: **%s**", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "---",
  "",
  "## Table of Contents",
  "",
  "1. [Final Genome Assemblies](#1-final-genome-assemblies)",
  "2. [Assembly QC Summary](#2-assembly-qc-summary)",
  "3. [Visual Summary](#3-visual-summary)",
  "4. [Assembly QC Across Pipeline Stages](#4-assembly-qc-across-pipeline-stages)",
  "5. [Telomere Detection](#5-telomere-detection)",
  "6. [Pairwise Alignment Summary](#6-pairwise-alignment-summary)",
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
  "Final gap-filled, chromosome-level assemblies produced by this pipeline.",
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
    sprintf("- **Interactive QC report (HTML):** [%s](%s)", e$filename, rel_path(e$subdir, e$filename)))
}
if (length(resource_links) > 0) {
  md <- c(md, "### Key Resources", "", resource_links, "")
}

# Assembly table
if (nrow(assemblies) > 0) {
  asm_table <- assemblies %>%
    arrange(id) %>%
    mutate(link = sprintf("[%s](%s)", filename, rel_path(subdir, filename))) %>%
    select(`Haplotype ID` = id, `Assembly File` = link)
  md <- c(md, "### Assembly Files", "", make_markdown_table(asm_table), "")
} else {
  md <- c(md, "*No final assemblies found in manifest.*", "")
}

# =============================================================================
# Section 2: Assembly QC Summary — per-haplotype, not averaged
# =============================================================================
md <- c(md, "## 2. Assembly QC Summary", "")

if (!is.null(qc_data) && nrow(qc_data) > 0) {

  # --- Per-haplotype overview table ---
  md <- c(md, "### Overview (Final Assemblies)", "")

  overview_long <- qc_data %>%
    filter(stage == "final") %>%
    pivot_longer(cols = c(hap1, hap2, both),
                 names_to = "haplotype",
                 values_to = "value",
                 values_drop_na = TRUE) %>%
    filter(metric %in% c("hifi_depth", "Total length", "L90", "auN",
                          "complete", "Largest contig", "GC (%)", "qv",
                          "kmer_completeness"))

  # For metrics that are diploid-level (qv, kmer_completeness), keep only "both"
  # For per-haplotype metrics, keep hap1 and hap2
  overview_filtered <- overview_long %>%
    filter(
      (metric %in% c("qv", "kmer_completeness") & haplotype == "both") |
      (!metric %in% c("qv", "kmer_completeness") & haplotype %in% c("hap1", "hap2"))
    )

  if (nrow(overview_filtered) > 0) {
    overview_fmt <- overview_filtered %>%
      mutate(value_fmt = case_when(
        analysis == "busco" ~ percent(value, accuracy = 0.1),
        metric == "hifi_depth" ~ comma(value, accuracy = 0.1),
        metric == "Total length" ~ comma(value, accuracy = 0.01, scale = 1/1e9, suffix = " Gb"),
        metric %in% c("Largest contig", "auN") ~ comma(value, accuracy = 0.1, scale = 1/1e6, suffix = " Mb"),
        metric %in% c("GC (%)", "kmer_completeness") ~ percent(value, scale = 1, accuracy = 0.1),
        metric == "qv" ~ comma(value, accuracy = 0.1),
        TRUE ~ comma(value, accuracy = 1)
      )) %>%
      select(sample_id, metric, haplotype, value_fmt) %>%
      pivot_wider(names_from = haplotype, values_from = value_fmt, values_fill = "") %>%
      rename(Sample = sample_id, Metric = metric) %>%
      rename_with(~ case_when(
        . == "complete"          ~ "BUSCO",
        . == "kmer_completeness" ~ "K-mer Comp.",
        . == "Largest contig"    ~ "Largest Scaffold",
        TRUE ~ .
      ))

    # Blank out repeated sample names
    overview_fmt <- overview_fmt %>%
      mutate(Sample = if_else(lag(Sample, default = "") == Sample, "", Sample))

    md <- c(md, make_markdown_table(overview_fmt), "")
  }

  # --- Detailed final assembly table (collapsible) ---
  final_detail <- qc_data %>%
    filter(stage == "final") %>%
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
      metric == "hifi_depth" ~ comma(value, accuracy = 0.1),
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
  "Each row shows one sample with snail plots, contact maps, and the within-sample",
  "haplotype-vs-haplotype dotplot.",
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

# Determine which columns to include
col_headers <- c("Sample")
if (has_snails) col_headers <- c(col_headers, "Hap1 Snail", "Hap2 Snail")
if (has_cmaps)  col_headers <- c(col_headers, "Hap1 Contact Map", "Hap2 Contact Map")
if (has_dots)   col_headers <- c(col_headers, "Hap1 vs Hap2 Dotplot")

if (has_snails || has_cmaps || has_dots) {

  img_w <- args$img_width

  # Start HTML table
  md <- c(md, "<table>", "<tr>")
  for (h in col_headers) {
    md <- c(md, sprintf("  <th>%s</th>", h))
  }
  md <- c(md, "</tr>")

  for (sid in all_sample_ids) {
    hap1_id <- paste0(sid, "_hap1")
    hap2_id <- paste0(sid, "_hap2")

    md <- c(md, "<tr>")

    # Sample name
    md <- c(md, sprintf('  <td><b>%s</b></td>', sid))

    # Snail plots
    if (has_snails) {
      for (hid in c(hap1_id, hap2_id)) {
        snail_row <- snail_plots %>% filter(id == hid)
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
      for (hid in c(hap1_id, hap2_id)) {
        hap_cmaps <- contact_maps %>% filter(id == hid)
        best <- pick_best_contact_map(hap_cmaps)
        if (!is.null(best) && nrow(best) > 0) {
          src <- rel_path(best$subdir[1], best$filename[1])
          md <- c(md, sprintf('  <td>%s</td>', img_tag(src, hid, width = img_w)))
        } else {
          md <- c(md, "  <td>—</td>")
        }
      }
    }

    # Within-sample dotplot (hap1 vs hap2)
    if (has_dots) {
      # Look for dotplot where id==hap1 and id2==hap2 (or vice versa)
      dot_row <- dotplots %>%
        filter(
          (id == hap1_id & id2 == hap2_id) |
          (id == hap2_id & id2 == hap1_id)
        )
      if (nrow(dot_row) > 0) {
        src <- rel_path(dot_row$subdir[1], dot_row$filename[1])
        md <- c(md, sprintf('  <td>%s</td>', img_tag(src, paste(hap1_id, "vs", hap2_id), width = img_w)))
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
        mutate(s1 = str_replace(id, "_hap[12]$", ""),
               s2 = str_replace(id2, "_hap[12]$", "")) %>%
        pull(s1) == dotplots %>%
        mutate(s2 = str_replace(id2, "_hap[12]$", "")) %>%
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

# =============================================================================
# Section 4: QC trend plots + cross-stage table
# =============================================================================
if (!is.null(qc_data) && nrow(qc_data) > 0 && nrow(qc_plots) > 0) {
  md <- c(md, "## 4. Assembly QC Across Pipeline Stages", "")

  plot_lines <- character()
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

  # Cross-stage comparison table
  cross_stage <- qc_data %>%
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
      metric == "hifi_depth" ~ comma(value, accuracy = 0.1),
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
}

# =============================================================================
# Section 5: Telomere Detection (collapsible)
# =============================================================================
telo_path <- args$telomere_summary
if (!str_detect(basename(telo_path), "NO_TELOMERES") &&
    file.exists(telo_path) && file.size(telo_path) > 0) {
  telo_data <- tryCatch(read_tsv(telo_path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(telo_data) && nrow(telo_data) > 0) {
    md <- c(md,
      "## 5. Telomere Detection", "",
      make_collapsible(
        make_markdown_table(telo_data %>% mutate(across(everything(), as.character))),
        "Click to expand: Telomere detection summary"
      )
    )
  }
}

# =============================================================================
# Section 6: Pairwise Alignment Summary (collapsible)
# =============================================================================
pw_path <- args$pairwise_summary
if (!str_detect(basename(pw_path), "NO_PAIRWISE") &&
    file.exists(pw_path) && file.size(pw_path) > 0) {
  pw_data <- tryCatch(read_tsv(pw_path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(pw_data) && nrow(pw_data) > 0) {
    md <- c(md,
      "## 6. Pairwise Alignment Summary", "",
      make_collapsible(
        make_markdown_table(pw_data %>% mutate(across(everything(), as.character))),
        "Click to expand: Pairwise alignment metrics"
      )
    )
  }
}

# =============================================================================
# Footer
# =============================================================================
md <- c(md, "---", "",
  sprintf("*Report generated on %s by the Genome Assembly Pipeline.*",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

writeLines(md, args$output)
message(sprintf("Report written to: %s (%d lines)", args$output, length(md)))