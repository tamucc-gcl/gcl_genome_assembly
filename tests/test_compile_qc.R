# Run: Rscript tests/test_compile_qc.R
# Exercises the production compiler body with synthetic stage/contact measurements.
suppressPackageStartupMessages(library(tidyverse))
source_path <- "r_scripts/compile_qc.R"
parse(file = source_path) |> invisible()
lines <- readLines(source_path, warn = FALSE)
start <- grep("^# Helper functions$", lines)[1]
compiler <- parse(text = paste(lines[start:length(lines)], collapse = "\n"))
root <- tempfile("qc-regression-")
dir.create(root)
stages <- c("contig", "contig_organelle_filtered", "contig_purged", "contig_corrected",
            "contig_decontam", "scaffold", "scaffold_corrected", "scaffold_decontam",
            "scaffold_round2", "chimera_broken", "gap_filled", "teloclip_extended", "final")
expected_stages <- c("ctg.base", "ctg.org", "ctg.purged", "ctg.cor", "ctg.deco",
                     "scaf.base", "scaf.cor", "scaf.deco", "scaf2", "chim.brk",
                     "gap_fill", "teloclip", "final")
for (name in c("assembly", "bam", "pairs", "out")) dir.create(file.path(root, name))
metrics <- tribble(
  ~analysis, ~metric, ~hap1, ~hap2, ~both,
  "quast", "# contigs", 10, 12, NA,
  "quast", "Total length", 10000000, 11000000, NA,
  "quast", "N50", 2000000, 2200000, NA,
  "quast", "L50", 2, 2, NA,
  "busco", "complete", 900, 910, NA,
  "busco", "single", 880, 890, NA,
  "busco", "duplicated", 20, 20, NA,
  "busco", "fragmented", 50, 40, NA,
  "busco", "missing", 50, 50, NA,
  "busco", "total_busco", 1000, 1000, NA,
  "merqury_qv", "qv", 60, 61, 60.5,
  "merqury_completeness", "kmer_completeness", 75, 76, 98,
  "mapped_reads", "coverage", 20, 21, NA)
for (stage in stages) {
  write_tsv(mutate(metrics, sample_id = "test", qc_label = stage, .before = everything()),
            file.path(root, "assembly", paste0(stage, ".tsv")))
}
contacts <- tribble(
  ~checkpoint, ~cis_pairs_contig, ~trans_pairs_contig, ~trans_to_cis_contig,
  ~cis_pairs_scaffold, ~trans_pairs_scaffold, ~trans_to_cis_scaffold,
  "contig_filtered", 100, 100, 1, NA, NA, NA,
  "scaffold_space", 100, 100, 1, 180, 20, 20/180,
  "scaffold_round2_filtered", 200, 16, .08, NA, NA, NA,
  "scaffold_round2_space", 200, 16, .08, 210, 6, 6/210)
contacts <- contacts %>% mutate(haplotype_id = "test_hap1", pairs_total = cis_pairs_contig + trans_pairs_contig,
                               parse_total_pairs = 1000, retention_pct = pairs_total/10)
for (i in seq_len(nrow(contacts))) write_tsv(contacts[i, ], file.path(root, "pairs", paste0(i, ".tsv")))
write_tsv(tibble(haplotype_id = "test_hap1", checkpoint = "scaffold_round2_raw_map",
                 bam_mapped_pct = 95), file.path(root, "bam", "map.tsv"))
env <- new.env()
env$args <- list(assembly_dir = file.path(root, "assembly"), bam_dir = file.path(root, "bam"),
                 pairs_dir = file.path(root, "pairs"), output_dir = file.path(root, "out"))
eval(compiler, envir = env)
result <- read_csv(file.path(root, "out", "assembly_qc_metrics.csv"), show_col_types = FALSE)
stopifnot(identical(unique(result$stage), expected_stages))
value <- function(stage, metric) result %>% filter(.data$stage == .env$stage, .data$metric == .env$metric) %>% pull(hap1)
stopifnot(isTRUE(all.equal(value("ctg.deco", "trans_to_cis"), 1)),
          isTRUE(all.equal(value("scaf.base", "trans_to_cis"), 20/180)),
          isTRUE(all.equal(value("scaf.deco", "trans_to_cis"), .08)),
          isTRUE(all.equal(value("scaf2", "trans_to_cis"), 6/210)),
          identical(value("scaf.deco", "bam_mapped_pct"), 95),
          length(value("final", "trans_to_cis")) == 0L)
# A repeated measurement must be rejected instead of silently taking the first file.
file.copy(file.path(root, "pairs", "1.tsv"), file.path(root, "pairs", "duplicate.tsv"))
failure <- tryCatch({eval(compiler, envir = env); NULL}, error = conditionMessage)
stopifnot(grepl("duplicate Hi-C haplotype/checkpoint", failure, fixed = TRUE))
unlink(root, recursive = TRUE)
cat("PASS: stage order, AGP contact counts, mapping coordinates, no invented final contacts and duplicate rejection\n")
