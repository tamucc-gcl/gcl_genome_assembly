if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  
  # Initialize variables
  input_dir <- NULL
  output_dir <- NULL
  sample_id <- NULL
  qc_label <- NULL
  
  # Parse arguments
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--input_dir") {
      input_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output_dir") {
      output_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--sample_id") {
      sample_id <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--qc_label") {
      qc_label <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
} else {
  input_dir <- 'C:/Users/jdsel/Downloads/summary/Sde-CBau_104_Ex1_qc_inputs'
  output_dir <- 'C:/Users/jdsel/Downloads/summary'
  sample_id <- 'Sde-CBau_104_Ex1'
  qc_label <- 'contig'
}

message("Input directory: ", input_dir)
message("Output directory: ", output_dir)
message("Sample ID: ", sample_id)
message("QC Label: ", qc_label, "\n")

#### Libraries ####
library(tidyverse)
library(jsonlite)

#### Helpers ####
# Normalize a haplotype label from an assembly/file name (Phase 2 haploid support).
#   hap1/hap2 -> hap1/hap2
#   primary   -> hap1   (a haploid/collapsed assembly lands in the hap1 column; hap2/both = NA)
#   both/Both -> both   (merqury combined line, diploid only)
norm_hap <- function(x) {
  h <- str_to_lower(str_extract(x, regex('hap[0-9]|primary|both', ignore_case = TRUE)))
  if_else(h == 'primary', 'hap1', h)
}

#### Debugging ####
list.dirs(input_dir)
list.files(input_dir, recursive = TRUE)

#### Get Data ####
quast_out <- file.path(input_dir, 'quast/transposed_report.tsv') %>%
  read_tsv(show_col_types = FALSE) %>%
  mutate(haplotype = norm_hap(Assembly),
         .before = everything(),
         .keep = 'unused') %>%
  pivot_longer(cols = -haplotype,
               names_to = 'metric') %>%
  pivot_wider(names_from = haplotype,
              values_from = value)

busco_out <- list.files(file.path(input_dir, 'busco'),
           pattern = 'short_summary.*json',
           recursive = TRUE, full.names = TRUE) %>%
  tibble(filepath  = .) %>%
  mutate(haplotype = norm_hap(filepath)) %>%
  rowwise(haplotype) %>%
  reframe(read_json(filepath) %>%
            as_tibble() %>%
            select(complete = C,
                   single = S, 
                   duplicated = D, 
                   fragmented = `F`, 
                   missing = M, 
                   total_busco = dataset_total_buscos)) %>%
  mutate(total_busco = as.integer(total_busco)) %>%
  pivot_longer(cols = -haplotype,
               names_to = 'metric') %>%
  pivot_wider(names_from = haplotype,
              values_from = value)

mapping_out <- list.files(file.path(input_dir, 'mapping'),
                          pattern = 'avg_depth.txt',
                          full.names = TRUE, recursive = TRUE) %>%
  tibble(file = .) %>%
  mutate(haplotype = norm_hap(file)) %>%
  rowwise(haplotype) %>%
  reframe(coverage = read_lines(file) %>% 
            as.numeric()) %>%
  pivot_longer(cols = -haplotype,
               names_to = 'metric') %>%
  pivot_wider(names_from = haplotype,
              values_from = value)
 
qv_out <- read_tsv(file.path(input_dir, 'merqury', str_c(sample_id, 'qv', sep = '.')),
                   col_names = c('assembly_id', 'error_kmer', 'total_assembly_kmer', 'qv', 'error_rate'),
                   show_col_types = FALSE) %>%
  mutate(haplotype = norm_hap(assembly_id),
         .before = everything(),
         .keep = 'unused') %>%
  pivot_longer(cols = -haplotype,
               names_to = 'metric') %>%
  pivot_wider(names_from = haplotype,
              values_from = value)

completeness_out <- read_tsv(file.path(input_dir, 'merqury', str_c(sample_id, 'completeness.stats', sep = '.')),
                         col_names = c('assembly_id', 'kmer_set', 'kmer_found', 'total_hifi_kmer', 'kmer_completeness'),
                         show_col_types = FALSE) %>%
  mutate(haplotype = norm_hap(assembly_id),
         .before = everything(),
         .keep = 'unused') %>%
  select(-kmer_set) %>%
  pivot_longer(cols = -haplotype,
               names_to = 'metric') %>%
  pivot_wider(names_from = haplotype,
              values_from = value)

#### Join it all together ####
joined_data <- bind_rows(.id = 'analysis',
                         mapped_reads  = mapping_out,
          quast = quast_out,
          busco = busco_out,
          merqury_qv = qv_out,
          merqury_completeness = completeness_out) %>%
  # Add sample_id and qc_label columns for downstream joining

  mutate(sample_id = sample_id,
         qc_label = qc_label,
         .before = everything())

# Guarantee a uniform hap1/hap2/both schema so downstream (compile_qc.R) is ploidy-agnostic.
# Haploid samples only produce a hap1 column here; hap2/both are added as NA.
for (col in c('hap1', 'hap2', 'both')) {
  if (!col %in% names(joined_data)) joined_data[[col]] <- NA
}
joined_data <- joined_data %>%
  select(sample_id, qc_label, analysis, metric, hap1, hap2, both)

# Output filename now includes qc_label to prevent collisions
output_filename <- str_c(sample_id, '_', qc_label, '_qc_summary.tsv')

write_tsv(joined_data,
          file.path(output_dir, output_filename))

message("Wrote output to: ", file.path(output_dir, output_filename))
