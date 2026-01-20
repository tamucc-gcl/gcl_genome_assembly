if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  
  # Initialize variables
  input_dir <- NULL
  output_dir <- NULL
  sample_id <- NULL
  
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
    } else {
      i <- i + 1
    }
  }
} else {
  input_dir <- 'C:/Users/jdsel/Downloads/summary/Sde-CBau_104_Ex1_qc_inputs'
  output_dir <- 'C:/Users/jdsel/Downloads/summary'
  sample_id <- 'Sde-CBau_104_Ex1'
}

message("Input directory:", input_dir)
message("Output directory:", output_dir)
message("Sample ID:", sample_id, "\n")

#### Libraries ####
library(tidyverse)
library(jsonlite)

#### Debugging ####
list.dirs(input_dir)
list.files(input_dir, recursive = TRUE)

#### Get Data ####
quast_out <- file.path(input_dir, 'quast/transposed_report.tsv') %>%
  read_tsv(show_col_types = FALSE) %>%
  mutate(haplotype = str_extract(Assembly, 'hap[0-9]'),
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
  mutate(haplotype = str_extract(filepath, 'hap[0-9]')) %>%
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
  mutate(haplotype = str_extract(file, 'hap[0-9]')) %>%
  rowwise(haplotype) %>%
  reframe(hifi_depth = read_lines(file) %>% 
            as.numeric()) %>%
  pivot_longer(cols = -haplotype,
               names_to = 'metric') %>%
  pivot_wider(names_from = haplotype,
              values_from = value)
 
qv_out <- read_tsv(file.path(input_dir, 'merqury', str_c(sample_id, 'qv', sep = '.')),
                   col_names = c('assembly_id', 'error_kmer', 'total_assembly_kmer', 'qv', 'error_rate'),
                   show_col_types = FALSE) %>%
  mutate(haplotype = str_extract(assembly_id, '(hap[0-9])|(Both)') %>%
           str_to_lower(),
         .before = everything(),
         .keep = 'unused') %>%
  pivot_longer(cols = -haplotype,
               names_to = 'metric') %>%
  pivot_wider(names_from = haplotype,
              values_from = value)

completeness_out <- read_tsv(file.path(input_dir, 'merqury', str_c(sample_id, 'completeness.stats', sep = '.')),
                         col_names = c('assembly_id', 'kmer_set', 'kmer_found', 'total_hifi_kmer', 'kmer_completeness'),
                         show_col_types = FALSE) %>%
  mutate(haplotype = str_extract(assembly_id, '(hap[0-9])|(both)'),
         .before = everything(),
         .keep = 'unused') %>%
  select(-kmer_set) %>%
  pivot_longer(cols = -haplotype,
               names_to = 'metric') %>%
  pivot_wider(names_from = haplotype,
              values_from = value)

#### Join it all together ####
joined_data <- bind_rows(.id = 'analysis',
          mapped_hifi = mapping_out,
          quast = quast_out,
          busco = busco_out,
          merqury_qv = qv_out,
          merqury_completeness = completeness_out)

write_tsv(joined_data,
          file.path(output_dir, 
                    str_c(sample_id, '_qc_summary.tsv')))
