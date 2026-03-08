#!/usr/bin/env Rscript

#' =============================================================================
#' RIPARIAN (RIBBON) PLOT FROM PAF FILE
#' =============================================================================
#' Generates a riparian / ribbon synteny plot from a pairwise genome alignment
#' PAF file using gggenomes.  Chromosomes are drawn as horizontal bars for each
#' haplotype (top = reference, bottom = query) and syntenic links are drawn as
#' ribbons coloured by the reference chromosome.
#'
#' Usage:
#'   Rscript riparian_paf.R \
#'       --paf  <paf_file>       \
#'       --ref_fai  <ref.fai>    \
#'       --query_fai <query.fai> \
#'       --ref  <ref_id>         \
#'       --query <query_id>      \
#'       --output <output.png>
#'
#' Dependencies: gggenomes (CRAN), ggplot2, dplyr, argparse
#' =============================================================================

suppressPackageStartupMessages({
  library(gggenomes)
  library(ggplot2)
  library(dplyr)
  library(argparse)
})

# =============================================================================
# CLI arguments
# =============================================================================
parser <- ArgumentParser(description = "Generate riparian plot from PAF alignment")
parser$add_argument("--paf",       required = TRUE,  help = "PAF file (can be gzipped)")
parser$add_argument("--ref_fai",   required = TRUE,  help = "FAI index for reference assembly")
parser$add_argument("--query_fai", required = TRUE,  help = "FAI index for query assembly")
parser$add_argument("--ref",       required = TRUE,  help = "Reference haplotype ID (top row)")
parser$add_argument("--query",     required = TRUE,  help = "Query haplotype ID (bottom row)")
parser$add_argument("--output",    required = TRUE,  help = "Output file (.png, .pdf, .svg)")
parser$add_argument("--width",     type = "double",  default = 14,    help = "Plot width (inches)")
parser$add_argument("--height",    type = "double",  default = 6,     help = "Plot height (inches)")
parser$add_argument("--dpi",       type = "integer", default = 300,   help = "PNG resolution")
parser$add_argument("--min_aln",   type = "integer", default = 50000, help = "Min alignment length to plot (bp)")
parser$add_argument("--min_seq",   type = "integer", default = 1000000, help = "Min sequence length to include (bp)")
parser$add_argument("--alpha",     type = "double",  default = 0.45,  help = "Ribbon transparency (0-1)")

args <- parser$parse_args()

# =============================================================================
# Read sequence lengths from FAI files
# =============================================================================
message("Reading FAI files...")

read_fai_as_seqs <- function(fai_path, bin_id) {
  fai <- read.table(fai_path, header = FALSE, sep = "\t",
                    col.names = c("seq_id", "length", "offset", "linebases", "linewidth"))
  fai %>%
    transmute(
      bin_id = bin_id,
      seq_id = seq_id,
      length = length
    )
}

ref_seqs   <- read_fai_as_seqs(args$ref_fai,   args$ref)
query_seqs <- read_fai_as_seqs(args$query_fai, args$query)

# Filter to sequences above minimum length
ref_seqs   <- ref_seqs   %>% filter(length >= args$min_seq)
query_seqs <- query_seqs %>% filter(length >= args$min_seq)

if (nrow(ref_seqs) == 0 || nrow(query_seqs) == 0) {
  message("WARNING: No sequences above min_seq threshold (", args$min_seq,
          " bp). Creating empty plot.")
  png(args$output, width = args$width, height = args$height,
      units = "in", res = args$dpi)
  plot.new()
  text(0.5, 0.5, paste("No sequences >=", args$min_seq, "bp"), cex = 1.5)
  dev.off()
  quit(save = "no", status = 0)
}

# Sort sequences by length descending for a cleaner layout
ref_seqs   <- ref_seqs   %>% arrange(desc(length))
query_seqs <- query_seqs %>% arrange(desc(length))

seqs <- bind_rows(ref_seqs, query_seqs)

message("  Reference sequences:  ", nrow(ref_seqs),
        " (", sum(ref_seqs$length) / 1e6, " Mb)")
message("  Query sequences:      ", nrow(query_seqs),
        " (", sum(query_seqs$length) / 1e6, " Mb)")

# =============================================================================
# Read PAF and convert to gggenomes links format
# =============================================================================
message("Reading PAF file: ", args$paf)

# gggenomes read_paf expects specific columns; read manually for robustness
paf_cols <- c("qname", "qlen", "qstart", "qend",
              "strand",
              "tname", "tlen", "tstart", "tend",
              "nmatch", "alen", "mapq")

paf_raw <- tryCatch({
  read.table(args$paf, header = FALSE, sep = "\t", fill = TRUE,
             comment.char = "", quote = "",
             col.names = c(paf_cols, paste0("extra", 1:20)))[, 1:12]
}, error = function(e) {
  # Try gzipped
  con <- gzfile(args$paf, "rt")
  on.exit(close(con))
  df <- read.table(con, header = FALSE, sep = "\t", fill = TRUE,
                   comment.char = "", quote = "",
                   col.names = c(paf_cols, paste0("extra", 1:20)))[, 1:12]
  df
})
names(paf_raw) <- paf_cols

message("  Total alignments read: ", nrow(paf_raw))

# Filter by alignment length
paf_filt <- paf_raw %>% filter(alen >= args$min_aln)
message("  Alignments after min_aln filter (>= ", args$min_aln, " bp): ", nrow(paf_filt))

# Keep only alignments to sequences we're plotting
kept_ref_seqs   <- ref_seqs$seq_id
kept_query_seqs <- query_seqs$seq_id

paf_filt <- paf_filt %>%
  filter(tname %in% kept_ref_seqs, qname %in% kept_query_seqs)

message("  Alignments on plotted sequences: ", nrow(paf_filt))

if (nrow(paf_filt) == 0) {
  message("WARNING: No alignments pass filters. Creating empty plot.")
  png(args$output, width = args$width, height = args$height,
      units = "in", res = args$dpi)
  plot.new()
  text(0.5, 0.5, "No alignments pass filters", cex = 1.5)
  dev.off()
  quit(save = "no", status = 0)
}

# Build gggenomes links table
# gggenomes expects: seq_id, start, end (genome 1) and seq_id2, start2, end2 (genome 2)
# For strand: if '-', swap start2/end2 so gggenomes draws twisted ribbons
links <- paf_filt %>%
  transmute(
    bin_id  = args$ref,
    seq_id  = tname,
    start   = tstart + 1L,   # PAF is 0-based; gggenomes is 1-based
    end     = tend,
    bin_id2 = args$query,
    seq_id2 = qname,
    start2  = ifelse(strand == "+", qstart + 1L, qend),
    end2    = ifelse(strand == "+", qend,         qstart + 1L),
    strand  = strand,
    alen    = alen,
    mapq    = mapq
  )

# =============================================================================
# Build gggenomes plot
# =============================================================================
message("Generating riparian plot...")

# Assign a colour palette based on reference chromosomes
n_ref_chr <- nrow(ref_seqs)

# Choose palette: use a qualitative palette that scales well
if (n_ref_chr <= 8) {
  chr_palette <- RColorBrewer::brewer.pal(max(3, n_ref_chr), "Set2")[1:n_ref_chr]
} else if (n_ref_chr <= 12) {
  chr_palette <- RColorBrewer::brewer.pal(n_ref_chr, "Set3")
} else {
  # For many chromosomes, generate from hcl
  chr_palette <- hcl.colors(n_ref_chr, palette = "Dark 3")
}
names(chr_palette) <- ref_seqs$seq_id

# Add ref chromosome colour to links
links$ref_chr <- links$seq_id

p <- gggenomes(seqs = seqs, links = links) +
  geom_seq(linewidth = 3, color = "grey70") +
  geom_bin_label(size = 4, fontface = "bold") +
  geom_seq_label(aes(label = seq_id), size = 2.2, vjust = -0.8, check_overlap = TRUE) +
  geom_link(aes(fill = ref_chr), alpha = args$alpha, linewidth = 0.1, color = "grey30") +
  scale_fill_manual(
    values = chr_palette,
    name   = paste(args$ref, "chromosome"),
    guide  = guide_legend(ncol = 2, override.aes = list(alpha = 0.8))
  ) +
  labs(
    title    = paste0("Synteny: ", args$ref, " vs ", args$query),
    subtitle = paste0("Alignments \u2265 ", format(args$min_aln, big.mark = ","), " bp  |  ",
                      "Sequences \u2265 ", format(args$min_seq, big.mark = ","), " bp"),
    caption  = "minimap2 asm5 alignment  |  gggenomes riparian plot"
  ) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    plot.caption  = element_text(size = 8, color = "grey50"),
    legend.position = "right"
  )

# =============================================================================
# Save
# =============================================================================
ext <- tools::file_ext(args$output)
if (ext == "png") {
  ggsave(args$output, p, width = args$width, height = args$height,
         dpi = args$dpi, bg = "white")
} else if (ext == "pdf") {
  ggsave(args$output, p, width = args$width, height = args$height, device = "pdf")
} else if (ext == "svg") {
  ggsave(args$output, p, width = args$width, height = args$height, device = "svg")
} else {
  ggsave(args$output, p, width = args$width, height = args$height,
         dpi = args$dpi, bg = "white")
}

message("Riparian plot saved to: ", args$output)