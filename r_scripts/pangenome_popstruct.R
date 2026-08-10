#!/usr/bin/env Rscript
# ======================================================================================
# pangenome_popstruct.R  (workstream D)
# Repo location: r_scripts/pangenome_popstruct.R
#
# Renders per-haplotype population-structure figures from an `odgi similarity` long-format
# pairwise table (grouped by PanSN haplotype):
#   - PCoA scatter (classical MDS of the distance matrix; the PCA analogue for distances)
#   - neighbour-joining tree (ape::nj)
# Each haplotype (incl. the reference) is one point / one leaf. Small-n tolerant: with < 3
# groups a labelled placeholder is written instead of failing. Always writes both PNGs.
#
# Usage:
#   Rscript pangenome_popstruct.R <similarity.tsv> <label> <outdir>
# ======================================================================================

suppressPackageStartupMessages({ library(ggplot2); library(ape) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: pangenome_popstruct.R <similarity.tsv> <label> <outdir>")
sim_f <- args[1]; label <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
op <- function(x) file.path(outdir, paste0(label, x))
theme_set(theme_bw(base_size = 12))

placeholder <- function(path, msg) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 5) +
    theme_void() + xlim(-1, 1) + ylim(-1, 1)
  ggsave(path, p, width = 7, height = 5, dpi = 150)
}

# ---- build the haplotype x haplotype distance matrix from the odgi similarity table ----
D <- NULL; groups <- character(0)
if (file.exists(sim_f) && file.size(sim_f) > 0) {
  s <- tryCatch(read.delim(sim_f, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (!is.null(s) && nrow(s) > 0 && ncol(s) >= 3) {
    nm <- names(s)
    # group-name columns: prefer group.a/group.b, else the first two columns
    ga_c <- if ("group.a" %in% nm) "group.a" else nm[1]
    gb_c <- if ("group.b" %in% nm) "group.b" else nm[2]
    # distance column: euclidean (from -d) preferred; else derive from a jaccard similarity
    if ("euclidean" %in% nm) {
      s$.d <- suppressWarnings(as.numeric(s[[ "euclidean" ]]))
    } else if (any(grepl("jaccard", nm, ignore.case = TRUE))) {
      jc <- nm[grepl("jaccard", nm, ignore.case = TRUE)][1]
      s$.d <- 1 - suppressWarnings(as.numeric(s[[jc]]))
    } else {
      # last resort: last numeric column, treated as a distance
      numc <- nm[vapply(s, is.numeric, logical(1))]
      s$.d <- if (length(numc)) suppressWarnings(as.numeric(s[[tail(numc, 1)]])) else NA_real_
    }
    groups <- sort(unique(c(as.character(s[[ga_c]]), as.character(s[[gb_c]]))))
    D <- matrix(0, length(groups), length(groups), dimnames = list(groups, groups))
    for (i in seq_len(nrow(s))) {
      a <- as.character(s[[ga_c]][i]); b <- as.character(s[[gb_c]][i]); d <- s$.d[i]
      if (!is.na(a) && !is.na(b) && !is.na(d)) { D[a, b] <- d; D[b, a] <- d }
    }
  }
}
# short leaf labels: drop a trailing PanSN hap suffix like '#0' only when it is uninformative
short <- function(x) x

# ---- PCoA scatter ----------------------------------------------------------------------
if (!is.null(D) && length(groups) >= 3) {
  co <- tryCatch(cmdscale(as.dist(D), k = 2, eig = TRUE), error = function(e) NULL)
  if (!is.null(co)) {
    pts <- as.data.frame(co$points); names(pts) <- c("Dim1", "Dim2")
    pts$hap <- short(rownames(pts))
    ev <- co$eig; ev[ev < 0] <- 0
    lab <- function(k) if (sum(ev) > 0) sprintf("PCoA %d (%.1f%%)", k, 100 * ev[k] / sum(ev)) else sprintf("PCoA %d", k)
    p <- ggplot(pts, aes(Dim1, Dim2, label = hap)) +
      geom_point(size = 3, colour = "#2c7fb8") + geom_text(vjust = -0.8, size = 3) +
      labs(title = paste0(label, " \u2014 haplotype ordination (PCoA)"), x = lab(1), y = lab(2))
    ggsave(op(".pca.png"), p, width = 7, height = 5, dpi = 150)
  } else placeholder(op(".pca.png"), "PCoA failed")
} else placeholder(op(".pca.png"), sprintf("ordination needs \u22653 haplotypes\n(have %d)", length(groups)))

# ---- NJ tree ---------------------------------------------------------------------------
if (!is.null(D) && length(groups) >= 3) {
  tr <- tryCatch(nj(as.dist(D)), error = function(e) NULL)
  if (!is.null(tr)) {
    tr$tip.label <- short(tr$tip.label)
    png(op(".njtree.png"), width = 7 * 150, height = 5 * 150, res = 150)
    plot(tr, type = "unrooted", lab4ut = "axial",
         main = paste0(label, " \u2014 haplotype NJ tree (graph similarity)"))
    dev.off()
  } else placeholder(op(".njtree.png"), "NJ tree failed")
} else placeholder(op(".njtree.png"), sprintf("NJ tree needs \u22653 haplotypes\n(have %d)", length(groups)))
