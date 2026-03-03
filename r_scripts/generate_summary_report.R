#!/usr/bin/env Rscript

#' =============================================================================
#' GENERATE PIPELINE SUMMARY REPORT
#' =============================================================================
#' Reads data passed natively via Nextflow channels (no directory scanning).
#' Produces a Markdown report with:
#'   - Links to final genome assemblies, compiled QC CSV, interactive HTML report
#'   - Assembly QC tables (final overview + full final detail)
#'   - QC trend plots from COMPILE_FINAL_QC (collapsible)
#'   - Cross-stage comparison table (collapsible within the plots section)
#'   - Snail plots, contact maps, dotplots with proper <img> sizing
#'   - Telomere and pairwise alignment summaries (collapsible)
#'
#' The report publishes to {outdir}/ directly, so all relative paths are
#' simply subdir/filename (e.g., snail_plots/sample_hap1_snail.svg).
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
parser$add_argument("--img_width",        default = 600, type = "integer", help = "Image display width in pixels")

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

#' Relative path from report at {outdir}/ to {outdir}/{subdir}/{filename}
#' Since the report is at the root, this is just subdir/filename
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
  "3. [Assembly QC Across Pipeline Stages](#3-assembly-qc-across-pipeline-stages)",
  "4. [Snail Plots](#4-snail-plots-final-assemblies)",
  "5. [Contact Maps](#5-contact-maps-final-assemblies)",
  "6. [Pairwise Dotplots](#6-pairwise-dotplots)",
  "7. [Telomere Detection](#7-telomere-detection)",
  "8. [Pairwise Alignment Summary](#8-pairwise-alignment-summary)",
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
# Section 2: Assembly QC Overview
# =============================================================================
md <- c(md, "## 2. Assembly QC Summary", "")

if (!is.null(qc_data) && nrow(qc_data) > 0) {

  # --- Compact overview table ---
  md <- c(md, "### Overview (Final Assemblies)", "")

  overview <- qc_data %>%
    filter(stage == "final") %>%
    pivot_longer(cols = c(hap1, hap2, both),
                 names_to = "haplotype",
                 values_to = "value",
                 values_drop_na = TRUE) %>%
    filter(metric %in% c("hifi_depth", "Total length", "L90", "auN",
                          "complete", "Largest contig", "GC (%)", "qv",
                          "kmer_completeness")) %>%
    mutate(value = case_when(
      metric %in% c("kmer_completeness", "qv") & haplotype != "both" ~ NA_real_,
      TRUE ~ value
    )) %>%
    summarise(.by = c(sample_id, metric),
              value = mean(value, na.rm = TRUE)) %>%
    pivot_wider(names_from = metric, values_from = value)

  if (nrow(overview) > 0) {
    overview <- overview %>%
      mutate(
        across(any_of(c("hifi_depth", "L90")),      ~ comma(., accuracy = 1)),
        across(any_of("Total length"),               ~ comma(., accuracy = 0.01, scale = 1/1e9, suffix = " Gb")),
        across(any_of(c("Largest contig", "auN")),   ~ comma(., accuracy = 0.1, scale = 1/1e6, suffix = " Mb")),
        across(any_of("complete"),                   ~ percent(., accuracy = 0.1)),
        across(any_of(c("GC (%)", "kmer_completeness")), ~ percent(., scale = 1, accuracy = 0.1)),
        across(any_of("qv"),                         ~ comma(., accuracy = 0.1))
      ) %>%
      rename(Sample = sample_id) %>%
      rename_with(~ case_when(
        . == "complete"          ~ "BUSCO",
        . == "kmer_completeness" ~ "K-mer Comp.",
        . == "Largest contig"    ~ "Largest Scaffold",
        TRUE ~ .
      ))

    md <- c(md, make_markdown_table(overview %>% mutate(across(everything(), as.character))), "")
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
        "Click to expand: Detailed final assembly metrics"
      )
    )
  }

  # =============================================================================
  # Section 3: QC trend plots + cross-stage table
  # =============================================================================
  md <- c(md, "## 3. Assembly QC Across Pipeline Stages", "")

  # QC trend plots from COMPILE_FINAL_QC in a collapsible section
  if (nrow(qc_plots) > 0) {
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

    # Cross-stage comparison table nested inside plots section
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
}

# =============================================================================
# Section 4: Snail Plots
# =============================================================================
if (nrow(snail_plots) > 0) {
  md <- c(md, "## 4. Snail Plots (Final Assemblies)", "")

  for (sid in all_sample_ids) {
    md <- c(md, sprintf("### %s", sid), "")

    sample_snails <- snail_plots %>%
      filter(str_detect(id, fixed(sid))) %>%
      arrange(id)

    if (nrow(sample_snails) >= 2) {
      md <- c(md, "<table><tr>",
        apply(sample_snails, 1, function(row) {
          src <- rel_path(row["subdir"], row["filename"])
          sprintf("<td>%s<br><em>%s</em></td>", img_tag(src, row["id"]), row["id"])
        }),
        "</tr></table>", "")
    } else if (nrow(sample_snails) == 1) {
      row <- sample_snails[1, ]
      src <- rel_path(row$subdir, row$filename)
      md <- c(md, img_tag(src, row$id), "", sprintf("*%s*", row$id), "")
    }
  }
}

# =============================================================================
# Section 5: Contact Maps
# =============================================================================
if (nrow(contact_maps) > 0) {
  md <- c(md, "## 5. Contact Maps (Final Assemblies)", "")

  for (sid in all_sample_ids) {
    md <- c(md, sprintf("### %s", sid), "")

    sample_cmaps <- contact_maps %>%
      filter(str_detect(id, fixed(sid))) %>%
      arrange(id, filename)

    # Group by haplotype so multiple resolutions are together
    for (hap_id in unique(sample_cmaps$id)) {
      hap_cmaps <- sample_cmaps %>% filter(id == hap_id)

      if (nrow(hap_cmaps) > 1) {
        # Multiple resolutions — show highest res prominently, rest collapsible
        md <- c(md, sprintf("#### %s", hap_id), "")

        # Show first (or pick a specific resolution) at full size
        first_row <- hap_cmaps[1, ]
        src <- rel_path(first_row$subdir, first_row$filename)
        md <- c(md, img_tag(src, hap_id), "")

        if (nrow(hap_cmaps) > 1) {
          other_lines <- character()
          for (j in seq_len(nrow(hap_cmaps))) {
            row <- hap_cmaps[j, ]
            src <- rel_path(row$subdir, row$filename)
            other_lines <- c(other_lines,
              sprintf("- [%s](%s)", row$filename, src))
          }
          md <- c(md,
            make_collapsible(other_lines,
                             sprintf("All contact map resolutions for %s", hap_id)))
        }
      } else {
        row <- hap_cmaps[1, ]
        src <- rel_path(row$subdir, row$filename)
        md <- c(md, sprintf("#### %s", hap_id), "",
          img_tag(src, hap_id), "")
      }
    }
  }
}

# =============================================================================
# Section 6: Pairwise Dotplots
# =============================================================================
if (nrow(dotplots) > 0) {
  md <- c(md, "## 6. Pairwise Dotplots", "")

  for (i in seq_len(nrow(dotplots))) {
    row <- dotplots[i, ]
    src <- rel_path(row$subdir, row$filename)
    md <- c(md,
      sprintf("### %s vs %s", row$id, row$id2), "",
      img_tag(src, sprintf("Dotplot: %s vs %s", row$id, row$id2)), "")
  }
}

# =============================================================================
# Section 7: Telomere Detection (collapsible)
# =============================================================================
telo_path <- args$telomere_summary
if (!str_detect(basename(telo_path), "NO_TELOMERES") &&
    file.exists(telo_path) && file.size(telo_path) > 0) {
  telo_data <- tryCatch(read_tsv(telo_path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(telo_data) && nrow(telo_data) > 0) {
    md <- c(md,
      "## 7. Telomere Detection", "",
      make_collapsible(
        make_markdown_table(telo_data %>% mutate(across(everything(), as.character))),
        "Click to expand: Telomere detection summary"
      )
    )
  }
}

# =============================================================================
# Section 8: Pairwise Alignment Summary (collapsible)
# =============================================================================
pw_path <- args$pairwise_summary
if (!str_detect(basename(pw_path), "NO_PAIRWISE") &&
    file.exists(pw_path) && file.size(pw_path) > 0) {
  pw_data <- tryCatch(read_tsv(pw_path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(pw_data) && nrow(pw_data) > 0) {
    md <- c(md,
      "## 8. Pairwise Alignment Summary", "",
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