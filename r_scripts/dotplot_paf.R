#!/usr/bin/env Rscript

#' =============================================================================
#' DOTPLOT FROM PAF FILE
#' =============================================================================
#' Generates a dotplot visualization from a pairwise genome alignment PAF file
#'
#' Usage:
#'   Rscript dotplot_paf.R --paf <paf_file> --ref <ref_id> --query <query_id> --output <output_file>
#'
#' Arguments:
#'   --paf     Path to PAF file (can be gzipped)
#'   --ref     Reference genome identifier (y-axis label)
#'   --query   Query genome identifier (x-axis label)
#'   --output  Output file path (supports .png, .pdf, .svg)
#'   --width   Plot width in inches (default: 10)
#'   --height  Plot height in inches (default: 10)
#'   --dpi     Resolution for PNG output (default: 300)
#'
#' Note: pafr is installed by the SETUP_PAFR process before any alignment jobs run
#' =============================================================================

suppressPackageStartupMessages({
  library(pafr)
  library(ggplot2)
  library(argparse)
})

# =============================================================================
# Parse command line arguments
# =============================================================================
parser <- ArgumentParser(description = "Generate dotplot from PAF alignment file")
parser$add_argument("--paf", required = TRUE, help = "Path to PAF file (can be gzipped)")
parser$add_argument("--ref", required = TRUE, help = "Reference genome ID (y-axis)")
parser$add_argument("--query", required = TRUE, help = "Query genome ID (x-axis)")
parser$add_argument("--output", required = TRUE, help = "Output file path (.png, .pdf, .svg)")
parser$add_argument("--width", type = "double", default = 10, help = "Plot width in inches")
parser$add_argument("--height", type = "double", default = 10, help = "Plot height in inches")
parser$add_argument("--dpi", type = "integer", default = 300, help = "Resolution for PNG output")
parser$add_argument("--min_align", type = "integer", default = 0, help = "Minimum alignment length to plot")

args <- parser$parse_args()

# =============================================================================
# Read PAF and generate dotplot
# =============================================================================
message("Reading PAF file: ", args$paf)

# An empty / record-less PAF is common and legitimate: cross-species or no-homology
# pairs produce zero alignments, and the upstream process writes a gzipped placeholder
# so the file always exists. pafr::read_paf() does NOT return a 0-row object for such
# input -- it aborts with "subscript out of bounds" while coercing columns that aren't
# there. So we must detect emptiness BEFORE calling read_paf(), and also trap any
# read_paf() error defensively, funneling both cases into the empty-plot branch below.
read_paf_lines <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

paf_lines  <- tryCatch(read_paf_lines(args$paf), error = function(e) character(0))
has_record <- any(nzchar(trimws(paf_lines)))

paf_data <- NULL
if (has_record) {
  paf_data <- tryCatch(
    read_paf(args$paf),
    error = function(e) {
      message("WARNING: read_paf() failed (", conditionMessage(e),
              "); treating as no alignments.")
      NULL
    }
  )
}

# Check if PAF has any alignments (empty file, placeholder, or unreadable)
if (is.null(paf_data) || nrow(paf_data) == 0) {
  message("WARNING: No alignments in PAF file. Creating empty plot.")

  # Create empty plot with message
  empty_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("No alignments between\n", args$query, "\nand\n", args$ref),
             size = 6, hjust = 0.5) +
    theme_void() +
    labs(title = paste(args$query, "vs", args$ref),
         subtitle = "No alignments passing filters")

  ggsave(args$output, empty_plot,
         width = args$width, height = args$height, dpi = args$dpi)

  message("Empty plot saved to: ", args$output)
  quit(status = 0)
}

# Filter by minimum alignment length if specified
if (args$min_align > 0) {
  paf_data <- subset(paf_data, alen >= args$min_align)
  message("Filtered to ", nrow(paf_data), " alignments >= ", args$min_align, " bp")
}

# A non-empty PAF can still be reduced to zero rows by the min_align filter -- guard
# again so the dotplot() call downstream never sees an empty object.
if (nrow(paf_data) == 0) {
  message("WARNING: No alignments remain after min_align filter. Creating empty plot.")

  empty_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("No alignments between\n", args$query, "\nand\n", args$ref),
             size = 6, hjust = 0.5) +
    theme_void() +
    labs(title = paste(args$query, "vs", args$ref),
         subtitle = paste0("No alignments >= ", args$min_align, " bp"))

  ggsave(args$output, empty_plot,
         width = args$width, height = args$height, dpi = args$dpi)

  message("Empty plot saved to: ", args$output)
  quit(status = 0)
}

message("Generating dotplot with ", nrow(paf_data), " alignments")

# Generate dotplot
pair_plot <- dotplot(paf_data,
                     label_seqs = FALSE,
                     dashes = FALSE,
                     xlab = args$query,
                     ylab = args$ref) +
  theme_classic() +
  theme(
    panel.background = element_rect(colour = 'black', fill = NA),
    panel.grid.major = element_line(colour = 'grey80', linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(title = paste(args$query, "vs", args$ref))

# =============================================================================
# Save plot
# =============================================================================
message("Saving plot to: ", args$output)

ggsave(args$output, pair_plot,
       width = args$width,
       height = args$height,
       dpi = args$dpi)

message("Dotplot saved successfully")
