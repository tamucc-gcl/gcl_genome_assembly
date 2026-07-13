#!/usr/bin/env Rscript

#' =============================================================================
#' COMPILE FINAL QC REPORT
#' =============================================================================
#' Aggregates all QC metrics from the genome assembly pipeline into a single
#' consolidated report.
#'
#' Inputs:
#'   --assembly_dir: Directory containing assembly QC summary TSVs
#'   --bam_dir: Directory containing Hi-C BAM metrics TSVs
#'   --pairs_dir: Directory containing Hi-C pairs metrics TSVs
#'   --output_dir: Directory for output files
#'
#' Outputs:
#'   - assembly_qc_metrics.csv: Consolidated QC metrics
#'   - *.png: QC trend plots
#' =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
})

# =============================================================================
# Parse command line arguments
# =============================================================================
parser <- ArgumentParser(description = "Compile final QC report from all pipeline stages")
parser$add_argument("--assembly_dir", required = TRUE,
                    help = "Directory containing assembly QC summary TSVs")
parser$add_argument("--bam_dir", required = TRUE,
                    help = "Directory containing Hi-C BAM metrics TSVs")
parser$add_argument("--pairs_dir", required = TRUE,
                    help = "Directory containing Hi-C pairs metrics TSVs")
parser$add_argument("--output_dir", required = TRUE,
                    help = "Output directory for compiled report")

args <- parser$parse_args()

# =============================================================================
# Helper functions
# =============================================================================

#' Read all TSV files from a directory and bind them together
read_all_tsvs <- function(dir_path, file_pattern = "*.tsv") {
  files <- list.files(dir_path, pattern = "\\.tsv$", full.names = TRUE)
  
  if (length(files) == 0) {
    message(sprintf("No TSV files found in %s", dir_path))
    return(NULL)
  }
  
  message(sprintf("Reading %d files from %s", length(files), dir_path))
  
  dfs <- lapply(files, function(f) {
    tryCatch({
      df <- read_tsv(f, show_col_types = FALSE)
      df$source_file <- basename(f)
      df
    }, error = function(e) {
      warning(sprintf("Failed to read %s: %s", f, e$message))
      NULL
    })
  })
  
  dfs <- dfs[!sapply(dfs, is.null)]
  
  if (length(dfs) == 0) {
    return(NULL)
  }
  
  bind_rows(dfs)
}

plot_dims <- function(plot, base_width = 6, base_height = 5) {
  build <- ggplot_build(plot)
  layout <- build$layout$layout
  
  ncol <- max(layout$COL)
  nrow <- max(layout$ROW)
  
  list(width = base_width * ncol, height = base_height * nrow)
}

# =============================================================================
# Read all input data
# =============================================================================
message("=== Reading Assembly QC Summaries ===")
assembly_qc <- read_all_tsvs(args$assembly_dir)

message("\n=== Reading BAM Metrics ===")
bam_metrics <- read_all_tsvs(args$bam_dir)

message("\n=== Reading Pairs Metrics ===")
pairs_metrics <- read_all_tsvs(args$pairs_dir)

# =============================================================================
# Process and combine data
# =============================================================================

message("\n=== Processing Data ===")

# Extract sample_id from haplotype_id for joining
if (!is.null(bam_metrics)) {
  bam_metrics <- bam_metrics %>%
    mutate(sample_id = str_replace(haplotype_id, "_(hap[12]|primary)$", ""),
           haplotype_id = str_extract(haplotype_id, 'hap[12]|primary'),
           haplotype_id = if_else(haplotype_id == 'primary', 'hap1', haplotype_id),
           .before = everything())
}

if (!is.null(pairs_metrics)) {
  pairs_metrics <- pairs_metrics %>%
    mutate(sample_id = str_replace(haplotype_id, "_(hap[12]|primary)$", ""),
           haplotype_id = str_extract(haplotype_id, 'hap[12]|primary'),
           haplotype_id = if_else(haplotype_id == 'primary', 'hap1', haplotype_id),
           .before = everything())
}


#### Join together for nice single output ####

# --- 1. Determine pipeline stages present ---
# --- Build stage levels and labels ---
# Finalized assembly is QC'd unconditionally as 'final' (ASSEMBLY_QC_FINAL); teloclip-extended
# is its own intermediate (ASSEMBLY_QC_TELOCLIP, pre-FINALIZE); gap_filled is 'gap_fill'.
# Absent stages (e.g. teloclip off, or short-read which has none of these) drop via fct_drop().
all_assembly_stages <- unique(assembly_qc$qc_label)
last_assembly_stage <- 'final'
message(sprintf("  Final assembly stage: %s", last_assembly_stage))

stage_levels <- c('contig', 'contig_mito_filtered', 'contig_purged',
                  'contig_corrected', 'contig_decontam',
                  'scaffold', 'scaffold_corrected', 'scaffold_round2',
                  'gap_filled', 'teloclip_extended', 'final')

stage_labels <- c('ctg.base', 'ctg.mito', 'ctg.purged', 'ctg.cor', 'ctg.deco',
                  'scaf.base', 'scaf.cor', 'scaf2',
                  'gap_fill', 'teloclip', 'final')

# --- 3. Process assembly QC ---
fixed_assembly <- assembly_qc %>%
  rename(stage = qc_label) %>%
  mutate(stage = factor(stage,
                        levels = stage_levels,
                        labels = stage_labels) %>%
           fct_drop(),
         analysis = case_when(str_detect(analysis, 'merqury') ~ 'merqury',
                              TRUE ~ analysis)) %>%
  filter(!is.na(stage)) %>%   # Drop any unrecognised qc_labels

  arrange(stage) %>%
  select(-source_file)

# --- Helper: resolve Hi-C checkpoint to assembly stage label ---
# The factor levels are now the LABELS (ctg.base, ctg.cor, scaf.base, etc.),
# so we must match against those, not the original qc_label names.
asm_stage_levels <- levels(fixed_assembly$stage)

last_ctg_stage  <- {
  idx <- str_which(asm_stage_levels, '^ctg\\.')
  if (length(idx) > 0) asm_stage_levels[max(idx)] else NA_character_
}
first_scaf_stage <- {
  idx <- str_which(asm_stage_levels, '^scaf')
  if (length(idx) > 0) asm_stage_levels[min(idx)] else NA_character_
}
last_scaf_stage <- {
  idx <- str_which(asm_stage_levels, '^scaf')
  if (length(idx) > 0) asm_stage_levels[max(idx)] else NA_character_
}
# scaf2 stage (for scaffold_round2_space/filtered): look for 'scaf2' explicitly,
# fall back to last scaffold stage
scaf2_stage <- {
  if ('scaf2' %in% asm_stage_levels) 'scaf2' else last_scaf_stage
}

message(sprintf("  Stage resolution: last_ctg=%s, first_scaf=%s, scaf2=%s, last_scaf=%s",
                last_ctg_stage, first_scaf_stage, scaf2_stage, last_scaf_stage))

# --- 4. Process BAM metrics ---
fixed_bam <- bam_metrics %>%
  select(-source_file) %>%
  pivot_longer(cols = where(is.numeric),
               names_to = 'metric') %>%
  pivot_wider(names_from = haplotype_id) %>%
  mutate(stage = case_when(
           checkpoint == 'contig_raw_map'          ~ last_ctg_stage,
           checkpoint == 'scaffold_round2_raw_map' ~ last_scaf_stage,
           checkpoint == 'final_raw_map'           ~ 'final'
         ) %>%
           factor(levels = asm_stage_levels),
         .keep = 'unused',
         .before = 'metric') %>%
  filter(!is.na(stage)) %>%   # Drop checkpoints that didn't resolve
  arrange(stage) %>%
  mutate(analysis = 'mapped_hic')

# --- 5. Process pairs metrics ---
fixed_pairs <- pairs_metrics %>%
  select(-source_file) %>%
  pivot_longer(cols = where(is.numeric),
               names_to = 'metric',
               values_drop_na = TRUE) %>%
  pivot_wider(names_from = haplotype_id) %>%
  filter(!metric %in% c('cis_pairs_scaffold', 'trans_pairs_scaffold', 'trans_to_cis_scaffold')) %>%
  mutate(metric = str_remove_all(metric, c('_contig|_scaffold')),
         metric = case_when(metric == 'parse_total_pairs' ~ 'mapped_pairs',
                            metric == 'pairs_total' ~ 'retained_pairs',
                            TRUE ~ metric),
         stage = case_when(
           checkpoint == 'contig_filtered'          ~ last_ctg_stage,
           checkpoint == 'scaffold_space'           ~ first_scaf_stage,
           checkpoint == 'scaffold_round2_space'    ~ scaf2_stage,
           checkpoint == 'scaffold_round2_filtered' ~ scaf2_stage,
           checkpoint == 'final_filtered'           ~ 'final'
         ) %>%
           factor(levels = asm_stage_levels),
         .keep = 'unused',
         .before = 'metric') %>%
  filter(!is.na(stage)) %>%   # Drop checkpoints that didn't resolve
  arrange(stage) %>%
  mutate(analysis = 'hic_contact')

# --- 6. Combine and write (deduplicate in case checkpoints overlap) ---
full_qc_data <- bind_rows(fixed_assembly, fixed_bam, fixed_pairs) %>%
  distinct(sample_id, stage, analysis, metric, .keep_all = TRUE) %>%
  arrange(stage, sample_id)

write_csv(full_qc_data,
          file.path(args$output_dir, 'assembly_qc_metrics.csv'))

#### Summary Plots ####
trans_cis_plot <- full_qc_data %>%
  filter(metric %in% c('trans_to_cis')) %>%
  pivot_longer(cols = c(hap1, hap2)) %>%
  ggplot(aes(x = stage, y = value,
             colour = sample_id,
             shape = name)) +
  geom_line(aes(group = interaction(sample_id, name))) +
  geom_point() +
  labs(shape = 'Haplotype',
       colour = 'Sample',
       y = 'HiC trans:cis ratio',
       x = 'Assembly Stage') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black', fill = NA),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.key = element_blank())

trans_cis_dims <- plot_dims(trans_cis_plot)
ggsave(file.path(args$output_dir,
                 'trans_cis.png'),
       plot = trans_cis_plot, 
       width = trans_cis_dims$width, 
       height = trans_cis_dims$height)

contigs_plot <- full_qc_data %>%
  filter(str_detect(metric, '# contigs')) %>%
  pivot_longer(cols = c(hap1, hap2)) %>%
  mutate(contig_size = str_extract(metric, '[0-9]+') %>% as.numeric() %>% replace_na(0),
         metric = str_replace(metric, '[0-9]+', scales::comma(contig_size)),
         metric = fct_reorder(metric, contig_size)) %>%
  ggplot(aes(x = stage, y = value,
             colour = sample_id,
             shape = name)) +
  geom_line(aes(group = interaction(sample_id, name))) +
  geom_point() +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~metric,
             scales = 'free_y',
             ncol = 2) +
  labs(shape = 'Haplotype',
       colour = 'Sample',
       y = 'Number of Contigs',
       x = 'Assembly Stage') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black', fill = NA),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.key = element_blank(),
        strip.background = element_blank())

contig_dims <- plot_dims(contigs_plot)
ggsave(file.path(args$output_dir,
                 'contig_count.png'),
       plot = contigs_plot, 
       width = contig_dims$width, 
       height = contig_dims$height)

size_plot <- full_qc_data %>%
  filter(str_detect(metric, 'Total length')) %>%
  pivot_longer(cols = c(hap1, hap2)) %>%
  mutate(contig_size = str_extract(metric, '[0-9]+') %>% as.numeric() %>% replace_na(0),
         metric = str_replace(metric, '[0-9]+', scales::comma(contig_size)),
         metric = fct_reorder(metric, contig_size)) %>%
  ggplot(aes(x = stage, y = value,
             colour = sample_id,
             shape = name)) +
  geom_line(aes(group = interaction(sample_id, name))) +
  geom_point() +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~metric,
             scales = 'free_y',
             ncol = 2) +
  labs(shape = 'Haplotype',
       colour = 'Sample',
       y = 'Assembly Size (bp)',
       x = 'Assembly Stage') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black', fill = NA),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.key = element_blank(),
        strip.background = element_blank())

size_dims <- plot_dims(size_plot)
ggsave(file.path(args$output_dir,
                 'contig_length.png'),
       plot = size_plot, 
       width = size_dims$width, 
       height = size_dims$height)

misc_quast_plots <- full_qc_data %>%
  filter(str_detect(metric, 'Largest contig|GC|N[0-9]0|L[0-9]0|auN|s per 100')) %>% 
  pivot_longer(cols = c(hap1, hap2)) %>%
  mutate(contig_size = str_extract(metric, '[0-9]+') %>% as.numeric() %>% replace_na(0),
         metric = str_replace(metric, '[0-9]+', scales::comma(contig_size)),
         metric = fct_reorder(metric, contig_size)) %>%
  ggplot(aes(x = stage, y = value,
             colour = sample_id,
             shape = name)) +
  geom_line(aes(group = interaction(sample_id, name))) +
  geom_point() +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~metric,
             scales = 'free_y',
             ncol = 2) +
  labs(shape = 'Haplotype',
       colour = 'Sample',
       y = 'value',
       x = 'Assembly Stage') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black', fill = NA),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.key = element_blank(),
        strip.background = element_blank())

misc_dims <- plot_dims(misc_quast_plots)
ggsave(file.path(args$output_dir,
                 'quast_misc.png'),
       plot = misc_quast_plots, 
       width = misc_dims$width, 
       height = misc_dims$height)

busco_plot <- full_qc_data %>%
  filter(analysis == 'busco') %>%
  pivot_longer(cols = c(hap1, hap2)) %>%
  mutate(value = value / value[metric == 'total_busco'],
         .by = c(sample_id, stage, name)) %>%
  filter(metric != 'total_busco') %>%
  mutate(metric = factor(metric, levels = c('complete',
                                            'single', 
                                            'duplicated',
                                            'fragmented', 
                                            'missing'))) %>%
  ggplot(aes(x = stage, y = value,
             colour = sample_id,
             shape = name)) +
  geom_line(aes(group = interaction(sample_id, name))) +
  geom_point() +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~metric,
             scales = 'free_y',
             ncol = 2) +
  labs(shape = 'Haplotype',
       colour = 'Sample',
       y = 'BUSCO Gene %',
       x = 'Assembly Stage') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black', fill = NA),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.key = element_blank(),
        strip.background = element_blank())

busco_dims <- plot_dims(busco_plot)
ggsave(file.path(args$output_dir,
                 'busco.png'),
       plot = busco_plot, 
       width = busco_dims$width, 
       height = busco_dims$height)

kmer_plot <- full_qc_data %>%
  filter(str_detect(metric, 'qv|kmer_completeness')) %>%
  pivot_longer(cols = c(hap1, hap2, both)) %>%
  ggplot(aes(x = stage, y = value,
             colour = sample_id,
             shape = name)) +
  geom_line(aes(group = interaction(sample_id, name))) +
  geom_point() +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~metric,
             scales = 'free_y',
             ncol = 2) +
  labs(shape = 'Haplotype',
       colour = 'Sample',
       y = 'value',
       x = 'Assembly Stage') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black', fill = NA),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.key = element_blank(),
        strip.background = element_blank())

kmer_dims <- plot_dims(kmer_plot)
ggsave(file.path(args$output_dir,
                 'kmer.png'),
       plot = kmer_plot, 
       width = kmer_dims$width, 
       height = kmer_dims$height)

message("\n=== QC Compilation Complete ===")