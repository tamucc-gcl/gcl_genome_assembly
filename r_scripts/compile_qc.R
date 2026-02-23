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
#'   - final_qc_report.tsv: Consolidated QC metrics
#'   - final_qc_report.html: Optional HTML visualization (if rmarkdown available)
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
# TODO: Implement your specific data processing logic here
#
# The assembly_qc dataframe will contain columns from COMBINE_ASSEMBLY_QC outputs
# The bam_metrics dataframe will contain columns:
#   - haplotype_id, checkpoint, bam_total_align, bam_mapped_align, 
#     bam_mapped_pct, bam_primary_align, bam_primary_mapped, bam_primary_mapped_pct
# The pairs_metrics dataframe will contain columns:
#   - haplotype_id, checkpoint, pairs_total, cis_pairs_contig, trans_pairs_contig,
#     trans_to_cis_contig, cis_pairs_scaffold, trans_pairs_scaffold, 
#     trans_to_cis_scaffold, parse_total_pairs, retention_pct

message("\n=== Processing Data ===")

# Example: Create a summary combining key metrics
# Modify this section based on your specific needs

# Extract sample_id from haplotype_id for joining
if (!is.null(bam_metrics)) {
  bam_metrics <- bam_metrics %>%
    mutate(sample_id = str_replace(haplotype_id, "_hap[12]$", ""),
           haplotype_id = str_extract(haplotype_id, 'hap[12]'),
           .before = everything())
}

if (!is.null(pairs_metrics)) {
  pairs_metrics <- pairs_metrics %>%
    mutate(sample_id = str_replace(haplotype_id, "_hap[12]$", ""),
           haplotype_id = str_extract(haplotype_id, 'hap[12]'),
           .before = everything())
}


#### Join together for nice single output ####
fixed_assembly <- assembly_qc %>%
  rename(stage = qc_label) %>%
  mutate(stage = factor(stage, 
                        levels = c('contig',
                                   'contig_purged',
                                   'contig_corrected',
                                   'contig_decontam',
                                   'scaffold',
                                   'scaffold_corrected',
                                   'scaffold_round2',
                                   'gap_filled')) %>%
           fct_drop(),
         analysis = case_when(str_detect(analysis, 'merqury') ~ 'merqury',
                              TRUE ~ analysis)) %>%
  arrange(stage) %>%
  select(-source_file)

fixed_bam <- bam_metrics %>%
  select(-source_file) %>%
  pivot_longer(cols = where(is.numeric),
               names_to = 'metric') %>%
  pivot_wider(names_from = haplotype_id) %>%
  mutate(stage = case_when(checkpoint == 'contig_raw_map' ~ levels(fixed_assembly$stage)[max(str_which(levels(fixed_assembly$stage), 'contig'))],
                           checkpoint == 'scaffold_round2_raw_map' ~ levels(fixed_assembly$stage)[max(str_which(levels(fixed_assembly$stage), 'scaffold'))],
                           checkpoint == 'final_raw_map' ~ 'gap_filled') %>%
           factor(levels = levels(fixed_assembly$stage)),
         .keep = 'unused',
         .before = 'metric') %>%
  arrange(stage) %>%
  mutate(analysis = 'mapped_hic')


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
         stage = case_when(checkpoint == 'contig_filtered' ~ levels(fixed_assembly$stage)[max(str_which(levels(fixed_assembly$stage), 'contig'))],
                           checkpoint == 'scaffold_space' ~ levels(fixed_assembly$stage)[min(str_which(levels(fixed_assembly$stage), 'scaffold'))],
                           checkpoint == 'scaffold_round2_space' ~ levels(fixed_assembly$stage)[min(str_which(levels(fixed_assembly$stage), 'scaffold')[-1])],
                           checkpoint == 'scaffold_round2_filtered' ~ levels(fixed_assembly$stage)[max(str_which(levels(fixed_assembly$stage), 'scaffold'))],
                           checkpoint == 'final_filtered' ~ 'gap_filled') %>%
           factor(levels = levels(fixed_assembly$stage)),
         .keep = 'unused',
         .before = 'metric') %>%
  arrange(stage) %>%
  mutate(analysis = 'hic_contact')

full_qc_data <- bind_rows(fixed_assembly,
          fixed_bam,
          fixed_pairs) %>%
  arrange(stage, sample_id) %>%
  mutate(stage = factor(stage, 
                        levels = c('contig',
                                   'contig_purged',
                                   'contig_corrected',
                                   'contig_decontam',
                                   'scaffold',
                                   'scaffold_corrected',
                                   'scaffold_round2',
                                   'gap_filled'),
                        labels = c('ctg.base',
                                   'ctg.purged',
                                   'ctg.cor',
                                   'ctg.deco',
                                   'scaf.base',
                                   'scaf.cor',
                                   'scaf2',
                                   'final')) %>%
           fct_drop())

write_csv(full_qc_data,
          file.path(args$output_dir,
                    'assembly_qc_metrics.csv'))

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
        # panel.grid.major = element_line(colour = 'grey80', linewidth = 0.25),
        # panel.grid.minor = element_blank(),
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
        # panel.grid.major = element_line(colour = 'grey80', linewidth = 0.25),
        # panel.grid.minor = element_blank(),
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
        # panel.grid.major = element_line(colour = 'grey80', linewidth = 0.25),
        # panel.grid.minor = element_blank(),
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
        # panel.grid.major = element_line(colour = 'grey80', linewidth = 0.25),
        # panel.grid.minor = element_blank(),
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
        # panel.grid.major = element_line(colour = 'grey80', linewidth = 0.25),
        # panel.grid.minor = element_blank(),
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
        # panel.grid.major = element_line(colour = 'grey80', linewidth = 0.25),
        # panel.grid.minor = element_blank(),
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