#!/usr/bin/env Rscript
# ======================================================================================
# pangenome_progressive.R  (workstream H)
# Repo location: r_scripts/pangenome_progressive.R
#
# Plots the empirical (incremental-construction) pangenome growth from PANGENOME_PROGRESSIVE:
# graph size (bp, and node count) as each assembly is added, reference-first.
#
# Usage: Rscript pangenome_progressive.R <progressive_growth.tsv> <label> <outdir>
# ======================================================================================

suppressPackageStartupMessages({ library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: pangenome_progressive.R <progressive_growth.tsv> <label> <outdir>")
tsv <- args[1]; label <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- file.path(outdir, paste0(label, ".progressive_growth.png"))
theme_set(theme_bw(base_size = 12))

d <- tryCatch(read.delim(tsv, header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
if (is.null(d) || nrow(d) < 1 || !all(c("k", "bp") %in% names(d))) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = "no progressive growth data", size = 5) + theme_void()
  ggsave(out, p, width = 7, height = 5, dpi = 150); quit(save = "no")
}
d$k <- as.integer(d$k); d$bp <- as.numeric(d$bp)

p <- ggplot(d, aes(k, bp)) +
  geom_line(linewidth = 0.9, colour = "#2c7fb8") +
  geom_point(size = 2.6, colour = "#2c7fb8") +
  scale_x_continuous(breaks = d$k) +
  scale_y_continuous(labels = function(x) paste0(round(x / 1e6), " Mb")) +
  labs(title = paste0(label, " \u2014 progressive pangenome growth (minigraph)"),
       subtitle = "graph size as each assembly is added (reference-first)",
       x = "assemblies (k)", y = "graph size")
ggsave(out, p, width = 7.5, height = 5, dpi = 150)
