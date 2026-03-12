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
        "5. [Mitochondrial Genome](#5-mitochondrial-genome)",
        "6. [Telomere Detection](#6-telomere-detection)",
        "7. [Pairwise Alignment Summary](#7-pairwise-alignment-summary)",
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
      sample_id = str_replace(id, "_hap[12]$", ""),
      hap = str_extract(id, "hap[12]$"),
      link = sprintf("[%s](%s)", filename, rel_path(subdir, filename))
    ) %>%
    select(sample_id, hap, link) %>%
    pivot_wider(names_from = hap, values_from = link, values_fill = "—") %>%
    rename(Sample = sample_id, `Haplotype 1` = hap1, `Haplotype 2` = hap2)
  
  md <- c(md, "### Assembly Files", "", make_markdown_table(asm_table), "")
} else {
  md <- c(md, "*No final assemblies found in manifest.*", "")
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
    filter(metric %in% c("hifi_depth", "Total length", "L90", "auN",
                         "complete", "Largest contig", "GC (%)", "qv",
                         "kmer_completeness"))
  
  if (nrow(overview_long) > 0) {
    # Format values
    fmt_value <- function(value, metric) {
      case_when(
        metric == "complete"                          ~ percent(value, accuracy = 0.1),
        metric == "hifi_depth"                        ~ comma(value, accuracy = 0.1),
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
      "hifi_depth"         = "HiFi Depth",
      "complete"           = "BUSCO Complete",
      "qv"                 = "QV",
      "kmer_completeness"  = "K-mer Completeness"
    )
    metric_order <- names(metric_display)
    
    # Diploid-level metrics (single value per sample)
    diploid_metrics <- c("qv", "kmer_completeness")
    
    # Per-haplotype metrics: format as "hap1 / hap2" slash notation
    per_hap <- overview_long %>%
      filter(!metric %in% diploid_metrics, haplotype %in% c("hap1", "hap2")) %>%
      mutate(value_fmt = fmt_value(value, metric)) %>%
      select(sample_id, metric, haplotype, value_fmt) %>%
      pivot_wider(names_from = haplotype, values_from = value_fmt, values_fill = "—") %>%
      mutate(combined = paste(hap1, "/", hap2)) %>%
      select(sample_id, metric, combined)
    
    # Diploid metrics: single value
    diploid <- overview_long %>%
      filter(metric %in% diploid_metrics, haplotype == "both") %>%
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
    
    md <- c(md, make_markdown_table(overview_wide), "")
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
    
    # Within-sample riparian plot (hap1 vs hap2)
    if (has_riparian) {
      rip_row <- riparian_plots %>%
        filter(
          (id == hap1_id & id2 == hap2_id) |
            (id == hap2_id & id2 == hap1_id)
        )
      if (nrow(rip_row) > 0) {
        src <- rel_path(rip_row$subdir[1], rip_row$filename[1])
        md <- c(md, sprintf('  <td>%s</td>', img_tag(src, paste(hap1_id, "vs", hap2_id, "riparian"), width = img_w)))
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
    
    # ---- FIX 1: Extract haplotype_id by stripping _scaffold_N suffix ----
    # Contig names follow the pattern {haplotype_id}_scaffold_{N},
    # e.g. "Sde-CBau_104_hap1_scaffold_2" or "Sde_CLim_110_hap2_scaffold_7".
    # The old regex (^[^_]+_[^_]+_hap[12]) assumed exactly two underscore-
    # delimited tokens before _hap, which fails when the haplotype_id itself
    # contains extra underscores (e.g. Sde_CLim_110_hap2).
    teloclip_summary <- teloclip_data %>%
      filter(extension_length > 0) %>%
      mutate(haplotype_id = str_remove(contig, "_scaffold_\\d+$")) %>%
      group_by(haplotype_id) %>%
      summarise(
        extensions = n(),
        total_bp_added = sum(extension_length),
        mean_extension_bp = round(mean(extension_length)),
        max_extension_bp = max(extension_length),
        .groups = "drop"
      )
    
    md <- c(md,
            "#### Telomere Extension (teloclip)", "",
            "Soft-clipped HiFi read overhangs containing telomeric motifs were used to",
            "extend scaffold ends missing telomeric sequence.", ""
    )
    
    if (nrow(teloclip_summary) > 0) {
      tc_table <- teloclip_summary %>%
        mutate(
          total_bp_added = scales::comma(total_bp_added),
          mean_extension_bp = scales::comma(mean_extension_bp),
          max_extension_bp = scales::comma(max_extension_bp)
        ) %>%
        rename(
          Haplotype = haplotype_id,
          Extensions = extensions,
          `Total bp Added` = total_bp_added,
          `Mean Extension (bp)` = mean_extension_bp,
          `Max Extension (bp)` = max_extension_bp
        )
      md <- c(md, make_markdown_table(tc_table), "")
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
# Footer
# =============================================================================
md <- c(md, "---", "",
        sprintf("*Report generated on %s by the Genome Assembly Pipeline.*",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

writeLines(md, args$output)
message(sprintf("Report written to: %s (%d lines)", args$output, length(md)))