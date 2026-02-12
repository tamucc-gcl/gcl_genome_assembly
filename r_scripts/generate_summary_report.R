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
message(sprintf("Looking in: %s", outdir))

# Get sample IDs from snail_plots directory (most reliable - these are final outputs)
snail_dir <- file.path(outdir, "snail_plots")
message(sprintf("Looking for snail plots in: %s", snail_dir))

if (dir.exists(snail_dir)) {
  all_snail_files <- list.files(snail_dir)
  message(sprintf("Snail dir contents (%d files):", length(all_snail_files)))
  message(paste("  ", all_snail_files, collapse = "\n"))
} else {
  message("WARNING: snail_dir does not exist!")
}

snail_files <- list.files(snail_dir, pattern = "_snail\\.svg$", full.names = FALSE)
message(sprintf("Found %d snail plot files", length(snail_files)))

# Extract haplotype IDs from snail plots, then get unique sample IDs
# Pattern: {sample_id}_hap{1,2}_{label}_snail.svg
haplotype_ids <- snail_files %>%
  str_extract("^.+_hap[12](?=_)")

sample_ids <- haplotype_ids %>%
  str_replace("_hap[12]$", "") %>%
  unique() %>%
  na.omit() %>%
  sort()

message(sprintf("Found %d samples: %s", length(sample_ids), paste(sample_ids, collapse = ", ")))

# If no snail plots, try to get samples from contact_maps or pairwise_alignments
if (length(sample_ids) == 0) {
  contact_dir <- file.path(outdir, "contact_maps")
  contact_files <- list.files(contact_dir, pattern = "_contact_map\\.png$", full.names = FALSE)
  haplotype_ids <- contact_files %>% str_extract("^.+_hap[12](?=_)")
  sample_ids <- haplotype_ids %>% str_replace("_hap[12]$", "") %>% unique() %>% na.omit() %>% sort()
  message(sprintf("From contact maps - found %d samples", length(sample_ids)))
}

samples <- sample_ids

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
  
  hap1_snail_link <- if (!is.na(row$hap1_snail)) sprintf("[snail](../snail_plots/%s)", row$hap1_snail) else "—"
  hap1_contact_link <- if (!is.na(row$hap1_contact)) sprintf("[contact](../contact_maps/%s)", row$hap1_contact) else "—"
  hap2_snail_link <- if (!is.na(row$hap2_snail)) sprintf("[snail](../snail_plots/%s)", row$hap2_snail) else "—"
  hap2_contact_link <- if (!is.na(row$hap2_contact)) sprintf("[contact](../contact_maps/%s)", row$hap2_contact) else "—"
  dotplot_link <- if (!is.na(row$dotplot)) sprintf("[dotplot](../pairwise_alignments/%s)", row$dotplot) else "—"
  
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s |",
                      row$sample_id, hap1_snail_link, hap1_contact_link, 
                      hap2_snail_link, hap2_contact_link, dotplot_link))
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
  
  hap1_snail_cell <- if (!is.na(row$hap1_snail)) sprintf('<img src="../snail_plots/%s">', row$hap1_snail) else "—"
  hap1_contact_cell <- if (!is.na(row$hap1_contact)) sprintf('<img src="../contact_maps/%s">', row$hap1_contact) else "—"
  hap2_snail_cell <- if (!is.na(row$hap2_snail)) sprintf('<img src="../snail_plots/%s">', row$hap2_snail) else "—"
  hap2_contact_cell <- if (!is.na(row$hap2_contact)) sprintf('<img src="../contact_maps/%s">', row$hap2_contact) else "—"
  dotplot_cell <- if (!is.na(row$dotplot)) sprintf('<img src="../pairwise_alignments/%s">', row$dotplot) else "—"
  
  html <- c(html, sprintf(
    '<tr><td><b>%s</b></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>',
    row$sample_id, hap1_snail_cell, hap1_contact_cell, hap2_snail_cell, hap2_contact_cell, dotplot_cell))
}

html <- c(html, "</table></body></html>")
writeLines(html, file.path(args$output_dir, "pipeline_summary_report.html"))

message("=== Done ===")