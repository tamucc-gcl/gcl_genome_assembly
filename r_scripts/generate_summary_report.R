#!/usr/bin/env Rscript

#' =============================================================================
#' GENERATE PIPELINE SUMMARY REPORT
#' =============================================================================
#' Creates comprehensive HTML and Markdown summary reports with:
#'   1. Visual table: snail plots (hap1, hap2) + within-sample dotplot per row
#'   2. QC metrics summary table (assembly, Hi-C mapping, pairs stats)
#'   3. Assembly QC plots (BUSCO, contiguity, QV)
#'   4. Links to all output files
#' =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
  library(knitr)
  library(kableExtra)
})

# =============================================================================
# Parse arguments
# =============================================================================
parser <- ArgumentParser(description = "Generate pipeline summary report")
parser$add_argument("--snail_dir", required = TRUE)
parser$add_argument("--dotplot_dir", required = TRUE)
parser$add_argument("--assembly_qc_dir", required = TRUE)
parser$add_argument("--output_dir", required = TRUE)
parser$add_argument("--outdir_base", default = "results")

args <- parser$parse_args()

# =============================================================================
# Helper functions
# =============================================================================
extract_sample_id <- function(haplotype_id) {
  str_replace(haplotype_id, "_hap[12]$", "")
}

extract_haplotype <- function(haplotype_id) {
  str_extract(haplotype_id, "hap[12]$")
}

img_to_base64 <- function(img_path) {
  if (!file.exists(img_path)) return("")
  ext <- tools::file_ext(img_path)
  mime <- switch(ext,
    "svg" = "image/svg+xml",
    "png" = "image/png",
    "jpg" = "image/jpeg",
    "jpeg" = "image/jpeg",
    "application/octet-stream"
  )
  raw <- readBin(img_path, "raw", file.info(img_path)$size)
  b64 <- base64enc::base64encode(raw)
  sprintf("data:%s;base64,%s", mime, b64)
}

# Alternative: use relative paths (for smaller HTML files)
img_relative_path <- function(img_path, base_dir = "summary_data") {
  if (!file.exists(img_path)) return("")
  file.path(base_dir, basename(dirname(img_path)), basename(img_path))
}

# =============================================================================
# 1. Load and organize snail plots
# =============================================================================
message("=== Loading snail plots ===")
snail_files <- list.files(args$snail_dir, pattern = "\\.(svg|png)$", full.names = TRUE)

snail_df <- tibble(path = snail_files) %>%
  mutate(
    filename = basename(path),
    # Expected format: {haplotype_id}_{qc_label}_snail.svg
    haplotype_id = str_extract(filename, "^[^_]+_[^_]+_hap[12]"),
    qc_label = str_extract(filename, "(?<=_hap[12]_)[^_]+(?=_snail)"),
    sample_id = extract_sample_id(haplotype_id),
    hap = extract_haplotype(haplotype_id)
  ) %>%
  filter(!is.na(haplotype_id))

message(sprintf("Found %d snail plots for %d samples", 
                nrow(snail_df), n_distinct(snail_df$sample_id)))

# =============================================================================
# 2. Load and organize dotplots
# =============================================================================
message("=== Loading dotplots ===")
dotplot_files <- list.files(args$dotplot_dir, pattern = "_dotplot\\.png$", full.names = TRUE)

dotplot_df <- tibble(path = dotplot_files) %>%
  mutate(
    filename = basename(path),
    # Expected format: {hap1_id}_vs_{hap2_id}_dotplot.png
    hap1_id = str_extract(filename, "^[^_]+_[^_]+_hap1"),
    hap2_id = str_extract(filename, "(?<=_vs_)[^_]+_[^_]+_hap2"),
    sample_id = extract_sample_id(hap1_id)
  ) %>%
  filter(!is.na(hap1_id))

message(sprintf("Found %d dotplots", nrow(dotplot_df)))

# =============================================================================
# 3. Load QC metrics
# =============================================================================
message("=== Loading QC metrics ===")
qc_files <- list.files(args$assembly_qc_dir, pattern = "\\.tsv$", full.names = TRUE)

# Read and combine all TSV files
all_qc <- map_dfr(qc_files, function(f) {
  tryCatch({
    df <- read_tsv(f, show_col_types = FALSE)
    df$source_file <- basename(f)
    df
  }, error = function(e) {
    warning(sprintf("Failed to read %s: %s", f, e$message))
    NULL
  })
})

# Parse assembly QC summaries (wide format with hap1, hap2, both columns)
assembly_qc <- all_qc %>%
  filter(str_detect(source_file, "qc_summary")) %>%
  select(-source_file) %>%
  distinct()

# Parse BAM metrics
bam_metrics <- all_qc %>%
  filter(str_detect(source_file, "bam_metrics|hic_bam")) %>%
  select(-source_file) %>%
  distinct()

# Parse pairs metrics
pairs_metrics <- all_qc %>%
  filter(str_detect(source_file, "pairs_metrics|hic_pairs")) %>%
  select(-source_file) %>%
  distinct()

message(sprintf("Loaded %d assembly QC records, %d BAM metrics, %d pairs metrics",
                nrow(assembly_qc), nrow(bam_metrics), nrow(pairs_metrics)))

# =============================================================================
# 4. Build visual comparison table (snail plots + dotplot)
# =============================================================================
message("=== Building visual comparison table ===")

# Get unique samples
samples <- sort(unique(c(snail_df$sample_id, dotplot_df$sample_id)))

# Choose embedding method based on file count
# Auto-switch to relative paths if >10 samples (avoids multi-MB HTML)
USE_BASE64 <- length(samples) <= 10

get_img_src <- function(img_path) {
  if (USE_BASE64) {
    img_to_base64(img_path)
  } else {
    # Copy to output and return relative path
    out_subdir <- file.path(args$output_dir, "images")
    dir.create(out_subdir, showWarnings = FALSE, recursive = TRUE)
    if (file.exists(img_path)) {
      file.copy(img_path, out_subdir, overwrite = TRUE)
      file.path("images", basename(img_path))
    } else {
      ""
    }
  }
}

visual_table_rows <- map_dfr(samples, function(sid) {
  # Get hap1 and hap2 snail plots for this sample
  hap1_snail <- snail_df %>% filter(sample_id == sid, hap == "hap1") %>% pull(path) %>% first()
  hap2_snail <- snail_df %>% filter(sample_id == sid, hap == "hap2") %>% pull(path) %>% first()
  
  # Get dotplot for this sample
  dotplot_path <- dotplot_df %>% filter(sample_id == sid) %>% pull(path) %>% first()
  
  tibble(
    sample_id = sid,
    hap1_snail_path = hap1_snail %||% NA_character_,
    hap2_snail_path = hap2_snail %||% NA_character_,
    dotplot_path = dotplot_path %||% NA_character_
  )
})

# =============================================================================
# 5. Build QC metrics summary table
# =============================================================================
message("=== Building QC metrics summary ===")

# Extract key metrics for final (gap_filled) stage
final_qc <- assembly_qc %>%
  filter(str_detect(stage, "gap_filled|final")) %>%
  select(sample_id, metric, hap1, hap2, both) %>%
  filter(metric %in% c(
    "total_length", "n50", "n_scaffolds", "largest_scaffold",
    "gc_percent", "qv", "kmer_completeness",
    "complete", "single", "duplicated", "fragmented", "missing"
  )) %>%
  pivot_wider(
    names_from = metric,
    values_from = c(hap1, hap2, both),
    names_glue = "{.value}_{metric}"
  )

# =============================================================================
# 6. Generate HTML report
# =============================================================================
message("=== Generating HTML report ===")

html_content <- c(
  '<!DOCTYPE html>',
  '<html lang="en">',
  '<head>',
  '  <meta charset="UTF-8">',
  '  <meta name="viewport" content="width=device-width, initial-scale=1.0">',
  '  <title>Genome Assembly Pipeline Summary Report</title>',
  '  <style>',
  '    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 40px; background: #f5f5f5; }',
  '    .container { max-width: 1600px; margin: 0 auto; background: white; padding: 40px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }',
  '    h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 15px; }',
  '    h2 { color: #34495e; margin-top: 40px; border-bottom: 2px solid #ecf0f1; padding-bottom: 10px; }',
  '    h3 { color: #7f8c8d; }',
  '    table { border-collapse: collapse; width: 100%; margin: 20px 0; }',
  '    th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }',
  '    th { background: #3498db; color: white; }',
  '    tr:nth-child(even) { background: #f9f9f9; }',
  '    tr:hover { background: #f1f1f1; }',
  '    .visual-table img { max-width: 300px; max-height: 300px; border: 1px solid #ddd; border-radius: 4px; }',
  '    .visual-table td { vertical-align: middle; text-align: center; }',
  '    .metric-good { color: #27ae60; font-weight: bold; }',
  '    .metric-warn { color: #f39c12; font-weight: bold; }',
  '    .metric-bad { color: #e74c3c; font-weight: bold; }',
  '    .links-section a { display: block; margin: 5px 0; color: #3498db; }',
  '    .links-section a:hover { color: #2980b9; }',
  '    .toc { background: #ecf0f1; padding: 20px; border-radius: 4px; margin-bottom: 30px; }',
  '    .toc ul { columns: 2; }',
  '    .toc a { color: #2c3e50; text-decoration: none; }',
  '    .toc a:hover { text-decoration: underline; }',
  '    .timestamp { color: #95a5a6; font-size: 0.9em; }',
  '    code { background: #ecf0f1; padding: 2px 6px; border-radius: 3px; }',
  '  </style>',
  '</head>',
  '<body>',
  '<div class="container">',
  sprintf('<h1>🧬 Genome Assembly Pipeline Summary Report</h1>'),
  sprintf('<p class="timestamp">Generated: %s</p>', format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  '',
  '<div class="toc">',
  '<h3>Table of Contents</h3>',
  '<ul>',
  '<li><a href="#visual-comparison">1. Visual Comparison (Snail Plots & Dotplots)</a></li>',
  '<li><a href="#qc-metrics">2. QC Metrics Summary</a></li>',
  '<li><a href="#assembly-stats">3. Assembly Statistics</a></li>',
  '<li><a href="#hic-stats">4. Hi-C Mapping Statistics</a></li>',
  '<li><a href="#output-links">5. Output File Links</a></li>',
  '</ul>',
  '</div>'
)

# Section 1: Visual comparison table
html_content <- c(html_content,
  '<h2 id="visual-comparison">1. Visual Comparison</h2>',
  '<p>Snail plots showing assembly contiguity and BUSCO completeness for each haplotype, ',
  'plus within-individual pairwise dotplots comparing hap1 vs hap2.</p>',
  '<table class="visual-table">',
  '<thead><tr><th>Sample</th><th>Haplotype 1 Snail Plot</th><th>Haplotype 2 Snail Plot</th><th>Hap1 vs Hap2 Dotplot</th></tr></thead>',
  '<tbody>'
)

for (i in seq_len(nrow(visual_table_rows))) {
  row <- visual_table_rows[i, ]
  
  hap1_img <- if (!is.na(row$hap1_snail_path) && file.exists(row$hap1_snail_path)) {
    sprintf('<img src="%s" alt="hap1 snail">', get_img_src(row$hap1_snail_path))
  } else { "<em>Not available</em>" }
  
  hap2_img <- if (!is.na(row$hap2_snail_path) && file.exists(row$hap2_snail_path)) {
    sprintf('<img src="%s" alt="hap2 snail">', get_img_src(row$hap2_snail_path))
  } else { "<em>Not available</em>" }
  
  dotplot_img <- if (!is.na(row$dotplot_path) && file.exists(row$dotplot_path)) {
    sprintf('<img src="%s" alt="dotplot">', get_img_src(row$dotplot_path))
  } else { "<em>Not available</em>" }
  
  html_content <- c(html_content,
    sprintf('<tr><td><strong>%s</strong></td><td>%s</td><td>%s</td><td>%s</td></tr>',
            row$sample_id, hap1_img, hap2_img, dotplot_img)
  )
}

html_content <- c(html_content, '</tbody></table>')

# Section 2: QC Metrics Summary
html_content <- c(html_content,
  '<h2 id="qc-metrics">2. QC Metrics Summary</h2>',
  '<p>Key quality metrics for final gap-filled assemblies.</p>'
)

if (nrow(final_qc) > 0) {
  # Build metrics table dynamically
  html_content <- c(html_content,
    '<table>',
    '<thead><tr><th>Sample</th><th>Haplotype</th><th>Total Length</th><th>N50</th><th>Scaffolds</th><th>QV</th><th>k-mer Comp.</th><th>BUSCO Complete</th></tr></thead>',
    '<tbody>'
  )
  
  for (i in seq_len(nrow(final_qc))) {
    row <- final_qc[i, ]
    
    # Hap1 row
    html_content <- c(html_content,
      sprintf('<tr><td rowspan="2"><strong>%s</strong></td><td>hap1</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>',
              row$sample_id,
              scales::comma(row$hap1_total_length %||% NA),
              scales::comma(row$hap1_n50 %||% NA),
              row$hap1_n_scaffolds %||% "—",
              sprintf("%.1f", row$hap1_qv %||% NA),
              sprintf("%.1f%%", (row$hap1_kmer_completeness %||% NA) * 100),
              sprintf("%.1f%%", (row$hap1_complete %||% NA) * 100))
    )
    
    # Hap2 row
    html_content <- c(html_content,
      sprintf('<tr><td>hap2</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>',
              scales::comma(row$hap2_total_length %||% NA),
              scales::comma(row$hap2_n50 %||% NA),
              row$hap2_n_scaffolds %||% "—",
              sprintf("%.1f", row$hap2_qv %||% NA),
              sprintf("%.1f%%", (row$hap2_kmer_completeness %||% NA) * 100),
              sprintf("%.1f%%", (row$hap2_complete %||% NA) * 100))
    )
  }
  
  html_content <- c(html_content, '</tbody></table>')
} else {
  html_content <- c(html_content, '<p><em>No final QC metrics available.</em></p>')
}

# Section 3: Assembly Statistics (full table)
html_content <- c(html_content,
  '<h2 id="assembly-stats">3. Detailed Assembly Statistics</h2>'
)

if (nrow(assembly_qc) > 0) {
  # Pivot to long format for display
  assembly_long <- assembly_qc %>%
    select(sample_id, stage, metric, hap1, hap2, both) %>%
    distinct()
  
  html_content <- c(html_content,
    '<details><summary>Click to expand full assembly statistics</summary>',
    '<table>',
    '<thead><tr><th>Sample</th><th>Stage</th><th>Metric</th><th>Hap1</th><th>Hap2</th><th>Both</th></tr></thead>',
    '<tbody>'
  )
  
  for (i in seq_len(nrow(assembly_long))) {
    row <- assembly_long[i, ]
    html_content <- c(html_content,
      sprintf('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>',
              row$sample_id, row$stage, row$metric,
              row$hap1 %||% "—", row$hap2 %||% "—", row$both %||% "—")
    )
  }
  
  html_content <- c(html_content, '</tbody></table></details>')
}

# Section 4: Hi-C Statistics
html_content <- c(html_content,
  '<h2 id="hic-stats">4. Hi-C Mapping Statistics</h2>'
)

if (nrow(pairs_metrics) > 0) {
  html_content <- c(html_content,
    '<details><summary>Click to expand Hi-C pairs statistics</summary>',
    '<table>',
    '<thead><tr>'
  )
  
  # Add column headers
  col_names <- names(pairs_metrics)
  for (col in col_names) {
    html_content <- c(html_content, sprintf('<th>%s</th>', col))
  }
  html_content <- c(html_content, '</tr></thead><tbody>')
  
  for (i in seq_len(nrow(pairs_metrics))) {
    html_content <- c(html_content, '<tr>')
    for (col in col_names) {
      val <- pairs_metrics[[col]][i]
      html_content <- c(html_content, sprintf('<td>%s</td>', val %||% "—"))
    }
    html_content <- c(html_content, '</tr>')
  }
  
  html_content <- c(html_content, '</tbody></table></details>')
} else {
  html_content <- c(html_content, '<p><em>No Hi-C pairs metrics available.</em></p>')
}

# Section 5: Output file links
html_content <- c(html_content,
  '<h2 id="output-links">5. Output File Links</h2>',
  '<div class="links-section">',
  '<h3>Final Assemblies</h3>',
  sprintf('<p>Location: <code>%s/assemblies/gap_filled/</code></p>', args$outdir_base),
  '',
  '<h3>QC Reports</h3>',
  '<ul>',
  sprintf('<li>Assembly QC: <code>%s/qc/assembly/</code></li>', args$outdir_base),
  sprintf('<li>Hi-C Mapping QC: <code>%s/qc/hic_mapping/</code></li>', args$outdir_base),
  sprintf('<li>QUAST Comparison: <code>%s/qc/quast_final/</code></li>', args$outdir_base),
  '</ul>',
  '',
  '<h3>Contact Maps</h3>',
  sprintf('<p>Location: <code>%s/contact_maps/</code></p>', args$outdir_base),
  '<p>View .mcool files with <a href="https://higlass.io/">HiGlass</a> or <a href="https://github.com/open2c/cooler">cooler</a></p>',
  '',
  '<h3>Pairwise Alignments</h3>',
  sprintf('<p>Location: <code>%s/pairwise_alignments/</code></p>', args$outdir_base),
  '',
  '<h3>Snail Plots</h3>',
  sprintf('<p>Location: <code>%s/snail_plots/</code></p>', args$outdir_base),
  '',
  '<h3>MultiQC Reports</h3>',
  sprintf('<p>Location: <code>%s/multiqc/</code></p>', args$outdir_base),
  '</div>'
)

# Close HTML
html_content <- c(html_content,
  '</div>',
  '</body>',
  '</html>'
)

# Write HTML
html_file <- file.path(args$output_dir, "pipeline_summary_report.html")
writeLines(html_content, html_file)
message(sprintf("Wrote HTML report: %s", html_file))

# =============================================================================
# 7. Generate Markdown report
# =============================================================================
message("=== Generating Markdown report ===")

md_content <- c(
  '# Genome Assembly Pipeline Summary Report',
  '',
  sprintf('**Generated:** %s', format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  '',
  '---',
  '',
  '## Table of Contents',
  '',
  '1. [Visual Comparison](#visual-comparison)',
  '2. [QC Metrics Summary](#qc-metrics-summary)',
  '3. [Assembly Statistics](#assembly-statistics)',
  '4. [Hi-C Mapping Statistics](#hi-c-mapping-statistics)',
  '5. [Output File Links](#output-file-links)',
  '',
  '---',
  '',
  '## Visual Comparison',
  '',
  'Snail plots and pairwise dotplots for each sample:',
  '',
  '| Sample | Hap1 Snail | Hap2 Snail | Hap1 vs Hap2 Dotplot |',
  '|--------|------------|------------|----------------------|'
)

for (i in seq_len(nrow(visual_table_rows))) {
  row <- visual_table_rows[i, ]
  hap1 <- if (!is.na(row$hap1_snail_path)) basename(row$hap1_snail_path) else "—"
  hap2 <- if (!is.na(row$hap2_snail_path)) basename(row$hap2_snail_path) else "—"
  dot <- if (!is.na(row$dotplot_path)) basename(row$dotplot_path) else "—"
  
  md_content <- c(md_content,
    sprintf('| %s | %s | %s | %s |', row$sample_id, hap1, hap2, dot)
  )
}

md_content <- c(md_content,
  '',
  '---',
  '',
  '## QC Metrics Summary',
  '',
  '### Final Assembly Quality',
  ''
)

if (nrow(final_qc) > 0) {
  md_content <- c(md_content,
    '| Sample | Hap | Total Length | N50 | Scaffolds | QV | k-mer Comp. | BUSCO |',
    '|--------|-----|-------------|-----|-----------|-----|-------------|-------|'
  )
  
  for (i in seq_len(nrow(final_qc))) {
    row <- final_qc[i, ]
    md_content <- c(md_content,
      sprintf('| %s | hap1 | %s | %s | %s | %.1f | %.1f%% | %.1f%% |',
              row$sample_id,
              scales::comma(row$hap1_total_length %||% NA),
              scales::comma(row$hap1_n50 %||% NA),
              row$hap1_n_scaffolds %||% "—",
              row$hap1_qv %||% NA,
              (row$hap1_kmer_completeness %||% NA) * 100,
              (row$hap1_complete %||% NA) * 100),
      sprintf('| | hap2 | %s | %s | %s | %.1f | %.1f%% | %.1f%% |',
              scales::comma(row$hap2_total_length %||% NA),
              scales::comma(row$hap2_n50 %||% NA),
              row$hap2_n_scaffolds %||% "—",
              row$hap2_qv %||% NA,
              (row$hap2_kmer_completeness %||% NA) * 100,
              (row$hap2_complete %||% NA) * 100)
    )
  }
}

md_content <- c(md_content,
  '',
  '---',
  '',
  '## Output File Links',
  '',
  sprintf('- **Final Assemblies:** `%s/assemblies/gap_filled/`', args$outdir_base),
  sprintf('- **Assembly QC:** `%s/qc/assembly/`', args$outdir_base),
  sprintf('- **Hi-C Mapping QC:** `%s/qc/hic_mapping/`', args$outdir_base),
  sprintf('- **Contact Maps:** `%s/contact_maps/`', args$outdir_base),
  sprintf('- **Pairwise Alignments:** `%s/pairwise_alignments/`', args$outdir_base),
  sprintf('- **Snail Plots:** `%s/snail_plots/`', args$outdir_base),
  sprintf('- **MultiQC:** `%s/multiqc/`', args$outdir_base)
)

# Write markdown
md_file <- file.path(args$output_dir, "pipeline_summary_report.md")
writeLines(md_content, md_file)
message(sprintf("Wrote Markdown report: %s", md_file))

message("\n=== Summary report generation complete ===")