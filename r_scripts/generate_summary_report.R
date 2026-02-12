#!/usr/bin/env Rscript

#' =============================================================================
#' GENERATE PIPELINE SUMMARY REPORT
#' =============================================================================
#' Scans published output directories and creates summary report
#' =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
})

parser <- ArgumentParser(description = "Generate pipeline summary report")
parser$add_argument("--outdir", required = TRUE, help = "Pipeline output directory")
parser$add_argument("--output_dir", default = ".", help = "Where to write reports")

args <- parser$parse_args()

outdir <- args$outdir

# =============================================================================
# Scan for files
# =============================================================================
message("=== Scanning output directories ===")

# Get sample IDs from QC TSVs
qc_dir <- file.path(outdir, "qc/assembly")
qc_files <- list.files(qc_dir, pattern = "_qc_summary\\.tsv$", 
                       full.names = TRUE, recursive = TRUE)

samples <- qc_files %>%
  basename() %>%
  str_extract("^[^_]+_[^_]+_[^_]+") %>%
  unique() %>%
  na.omit() %>%
  sort()

message(sprintf("Found %d samples: %s", length(samples), paste(samples, collapse = ", ")))

# Build expected filenames
visual_data <- tibble(sample_id = samples) %>%
  mutate(
    hap1_snail = sprintf("%s_hap1_gap_filled_snail.svg", sample_id),
    hap2_snail = sprintf("%s_hap2_gap_filled_snail.svg", sample_id),
    hap1_contact = sprintf("%s_hap1_final_1000000bp_contact_map.png", sample_id),
    hap2_contact = sprintf("%s_hap2_final_1000000bp_contact_map.png", sample_id),
    dotplot = sprintf("%s_hap1_vs_%s_hap2_dotplot.png", sample_id, sample_id)
  )

# =============================================================================
# Generate Markdown
# =============================================================================
message("=== Generating report ===")

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

md <- c(md,
  "", "---", "",
  "## 2. Output Locations",
  "",
  sprintf("- **Final Assemblies:** `%s/assemblies/gap_filled/`", outdir),
  sprintf("- **Snail Plots:** `%s/snail_plots/`", outdir),
  sprintf("- **Contact Maps:** `%s/contact_maps/`", outdir),
  sprintf("- **Pairwise Alignments:** `%s/pairwise_alignments/`", outdir),
  sprintf("- **Assembly QC:** `%s/qc/assembly/`", outdir)
)

writeLines(md, file.path(args$output_dir, "pipeline_summary_report.md"))

# =============================================================================
# Generate HTML with embedded images
# =============================================================================
html <- c(
  "<!DOCTYPE html>",
  "<html><head>",
  "<title>Pipeline Summary</title>",
  "<style>",
  "body{font-family:sans-serif;max-width:1600px;margin:40px auto;padding:20px}",
  "table{border-collapse:collapse;width:100%}",
  "th,td{border:1px solid #ddd;padding:8px;text-align:center}",
  "th{background:#4a90d9;color:white}",
  "img{max-width:250px;max-height:250px}",
  "</style>",
  "</head><body>",
  "<h1>Genome Assembly Pipeline Summary</h1>",
  sprintf("<p><em>Generated: %s</em></p>", format(Sys.time())),
  "<h2>Visual Comparison</h2>",
  "<table>",
  "<tr><th>Sample</th><th>Hap1 Snail</th><th>Hap1 Contact</th><th>Hap2 Snail</th><th>Hap2 Contact</th><th>Dotplot</th></tr>"
)

for (i in seq_len(nrow(visual_data))) {
  row <- visual_data[i, ]
  html <- c(html, sprintf(
    '<tr><td><b>%s</b></td><td><img src="../snail_plots/%s"></td><td><img src="../contact_maps/%s"></td><td><img src="../snail_plots/%s"></td><td><img src="../contact_maps/%s"></td><td><img src="../pairwise_alignments/%s"></td></tr>',
    row$sample_id, row$hap1_snail, row$hap1_contact, row$hap2_snail, row$hap2_contact, row$dotplot))
}

html <- c(html, "</table></body></html>")
writeLines(html, file.path(args$output_dir, "pipeline_summary_report.html"))

message("=== Done ===")