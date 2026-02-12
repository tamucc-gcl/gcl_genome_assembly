#!/usr/bin/env Rscript

#' =============================================================================
#' GENERATE PIPELINE SUMMARY REPORT
#' =============================================================================
#' Creates markdown and HTML reports with:
#'   1. Visual table with links to snail plots, contact maps, dotplots
#'   2. QC metrics table from compiled TSVs
#'   3. Links to output directories
#'
#' Scans published output directories for images
#' =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
})

parser <- ArgumentParser(description = "Generate pipeline summary report")
parser$add_argument("--qc_dir", required = TRUE, help = "Directory with compiled QC TSVs")
parser$add_argument("--outdir_base", required = TRUE, help = "Pipeline output base directory")
parser$add_argument("--output_dir", default = ".", help = "Where to write reports")

args <- parser$parse_args()

# =============================================================================
# Helper functions
# =============================================================================
extract_sample_id <- function(haplotype_id) {

  str_replace(haplotype_id, "_hap[12]$", "")
}

extract_hap <- function(haplotype_id) {
  str_extract(haplotype_id, "hap[12]$")
}

# =============================================================================
# Load QC data and derive sample/haplotype IDs
# =============================================================================
message("=== Loading QC data ===")

outdir <- args$outdir_base

assembly_qc_dir <- file.path(args$qc_dir, "assembly")
assembly_files <- list.files(assembly_qc_dir, pattern = "\\.tsv$", full.names = TRUE)

assembly_qc <- map_dfr(assembly_files, function(f) {
  tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
})

# Get unique sample IDs from QC data
if (nrow(assembly_qc) > 0 && "sample_id" %in% names(assembly_qc)) {
  samples <- sort(unique(assembly_qc$sample_id))
  message(sprintf("Found %d samples in QC data: %s", length(samples), paste(samples, collapse = ", ")))
} else {
  # Try to extract from filenames
  samples <- assembly_files %>%
    basename() %>%
    str_extract("^[^_]+_[^_]+_[^_]+") %>%
    unique() %>%
    na.omit()
  message(sprintf("Extracted %d samples from filenames", length(samples)))
}

# Build expected filenames based on sample IDs
# These are relative paths from the reports/ directory
visual_data <- tibble(sample_id = samples) %>%
  mutate(
    hap1_snail = sprintf("%s_hap1_gap_filled_snail.svg", sample_id),
    hap2_snail = sprintf("%s_hap2_gap_filled_snail.svg", sample_id),
    hap1_contact = sprintf("%s_hap1_final_1000000bp_contact_map.png", sample_id),
    hap2_contact = sprintf("%s_hap2_final_1000000bp_contact_map.png", sample_id),
    dotplot = sprintf("%s_hap1_vs_%s_hap2_dotplot.png", sample_id, sample_id)
  )

message(sprintf("Generated links for %d samples", nrow(visual_data)))

# =============================================================================
# Generate Markdown Report
# =============================================================================
message("=== Generating Markdown report ===")

md <- c(
  "# Genome Assembly Pipeline Summary Report",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "---",
  "",
  "## 1. Visual Comparison",
  "",
  "| Sample | Hap1 Snail | Hap1 Contact | Hap2 Snail | Hap2 Contact | Hap1 vs Hap2 |",
  "|--------|------------|--------------|------------|--------------|--------------|"
)

for (i in seq_len(nrow(visual_data))) {
  row <- visual_data[i, ]
  
  md <- c(md, sprintf("| %s | [snail](../snail_plots/%s) | [contact](../contact_maps/%s) | [snail](../snail_plots/%s) | [contact](../contact_maps/%s) | [dotplot](../pairwise_alignments/%s) |",
                      row$sample_id, row$hap1_snail, row$hap1_contact, 
                      row$hap2_snail, row$hap2_contact, row$dotplot))
}

# QC Metrics table
md <- c(md, "", "---", "", "## 2. QC Metrics Summary", "")

if (nrow(assembly_qc) > 0 && "stage" %in% names(assembly_qc)) {
  # Get final stage metrics
  final_qc <- assembly_qc %>%
    filter(str_detect(stage, "gap_filled|final")) %>%
    filter(metric %in% c("total_length", "n50", "n_scaffolds", "qv", "kmer_completeness", "complete"))
  
  if (nrow(final_qc) > 0) {
    md <- c(md,
      "| Sample | Metric | Hap1 | Hap2 |",
      "|--------|--------|------|------|"
    )
    
    for (i in seq_len(nrow(final_qc))) {
      row <- final_qc[i, ]
      md <- c(md, sprintf("| %s | %s | %s | %s |",
                          row$sample_id, row$metric,
                          row$hap1 %||% "—", row$hap2 %||% "—"))
    }
  }
} else {
  md <- c(md, "*No QC metrics available*")
}

# Output links
md <- c(md,
  "", "---", "",
  "## 3. Output File Locations",
  "",
  sprintf("- **Final Assemblies:** `%s/assemblies/gap_filled/`", outdir),
  sprintf("- **Snail Plots:** `%s/snail_plots/`", outdir),
  sprintf("- **Contact Maps:** `%s/contact_maps/`", outdir),
  sprintf("- **Pairwise Alignments:** `%s/pairwise_alignments/`", outdir),
  sprintf("- **Assembly QC:** `%s/qc/assembly/`", outdir),
  sprintf("- **Hi-C QC:** `%s/qc/hic_mapping/`", outdir),
  sprintf("- **QUAST Final:** `%s/qc/assembly/quast_final/`", outdir)
)

writeLines(md, file.path(args$output_dir, "pipeline_summary_report.md"))
message("Wrote pipeline_summary_report.md")

# =============================================================================
# Generate simple HTML (just wrap markdown table in basic HTML)
# =============================================================================
html <- c(
  "<!DOCTYPE html>",
  "<html><head>",
  "<title>Pipeline Summary Report</title>",
  "<style>",
  "body { font-family: sans-serif; max-width: 1400px; margin: 40px auto; padding: 20px; }",
  "table { border-collapse: collapse; width: 100%; margin: 20px 0; }",
  "th, td { border: 1px solid #ddd; padding: 10px; text-align: left; }",
  "th { background: #4a90d9; color: white; }",
  "tr:nth-child(even) { background: #f9f9f9; }",
  "a { color: #4a90d9; }",
  "img { max-width: 250px; max-height: 250px; }",
  "</style>",
  "</head><body>",
  "<h1>Genome Assembly Pipeline Summary Report</h1>",
  sprintf("<p><em>Generated: %s</em></p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "<hr>",
  "<h2>1. Visual Comparison</h2>",
  "<table>",
  "<tr><th>Sample</th><th>Hap1 Snail</th><th>Hap1 Contact</th><th>Hap2 Snail</th><th>Hap2 Contact</th><th>Hap1 vs Hap2</th></tr>"
)

for (i in seq_len(nrow(visual_data))) {
  row <- visual_data[i, ]
  
  html <- c(html, sprintf('<tr><td><strong>%s</strong></td><td><a href="../snail_plots/%s"><img src="../snail_plots/%s"></a></td><td><a href="../contact_maps/%s"><img src="../contact_maps/%s"></a></td><td><a href="../snail_plots/%s"><img src="../snail_plots/%s"></a></td><td><a href="../contact_maps/%s"><img src="../contact_maps/%s"></a></td><td><a href="../pairwise_alignments/%s"><img src="../pairwise_alignments/%s"></a></td></tr>',
                          row$sample_id, 
                          row$hap1_snail, row$hap1_snail,
                          row$hap1_contact, row$hap1_contact,
                          row$hap2_snail, row$hap2_snail,
                          row$hap2_contact, row$hap2_contact,
                          row$dotplot, row$dotplot))
}

html <- c(html, "</table>", "</body></html>")

writeLines(html, file.path(args$output_dir, "pipeline_summary_report.html"))
message("Wrote pipeline_summary_report.html")

message("\n=== Done ===")