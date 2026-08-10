#!/usr/bin/env Rscript
# ======================================================================================
# pangenome_popstruct.R  (workstream D)
# Repo location: r_scripts/pangenome_popstruct.R
#
# From an `odgi similarity -D '#' -p 2 -d` long-format table (grouped by PanSN haplotype),
# renders per-HAPLOTYPE and per-INDIVIDUAL ordination (PCoA) + neighbour-joining trees.
# Distances come from shared graph content (SNPs + indels + SVs). Small-n tolerant: a level
# with < 3 units gets a labelled placeholder. Always writes all four PNGs.
#
# Usage: Rscript pangenome_popstruct.R <similarity.tsv> <label> <outdir>
# ======================================================================================

suppressPackageStartupMessages({ library(ggplot2); library(ape); library(ggrepel) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: pangenome_popstruct.R <similarity.tsv> <label> <outdir>")
sim_f <- args[1]; label <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
op  <- function(x) file.path(outdir, paste0(label, x))
theme_set(theme_bw(base_size = 12))
DASH <- "\u2014"; GE <- "\u2265"

# short, readable labels: drop the common 'Sde-' prefix, and collapse the reference-style
# '_dip_hapN#0' suffix to '#N' so every haplotype reads uniformly (e.g. CBau_104#1, CMat_203#2)
short <- function(x) sub("_dip_hap([0-9]+)#0$", "#\\1", sub("^Sde-", "", x))

placeholder <- function(path, msg) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 5) +
    theme_void() + xlim(-1, 1) + ylim(-1, 1)
  ggsave(path, p, width = 7.5, height = 5.5, dpi = 150)
}

# ---- build the haplotype x haplotype distance matrix -----------------------------------
D <- NULL; groups <- character(0)
if (file.exists(sim_f) && file.size(sim_f) > 0) {
  s <- tryCatch(read.delim(sim_f, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (!is.null(s) && nrow(s) > 0 && ncol(s) >= 3) {
    nm <- names(s)
    ga_c <- if ("group.a" %in% nm) "group.a" else nm[1]
    gb_c <- if ("group.b" %in% nm) "group.b" else nm[2]
    dcands <- nm[grepl("distance$", nm, ignore.case = TRUE) & !grepl("^group", nm, ignore.case = TRUE)]
    dcol <- if ("jaccard.distance" %in% nm) "jaccard.distance"
            else if ("euclidean.distance" %in% nm) "euclidean.distance"
            else if (length(dcands)) dcands[1] else NA
    if (!is.na(dcol)) {
      s$.d <- suppressWarnings(as.numeric(s[[dcol]]))
      groups <- sort(unique(c(as.character(s[[ga_c]]), as.character(s[[gb_c]]))))
      D <- matrix(0, length(groups), length(groups), dimnames = list(groups, groups))
      for (i in seq_len(nrow(s))) {
        a <- as.character(s[[ga_c]][i]); b <- as.character(s[[gb_c]][i]); d <- s$.d[i]
        if (!is.na(a) && !is.na(b) && !is.na(d)) { D[a, b] <- d; D[b, a] <- d }
      }
    }
  }
}

# ---- plot one level (PCoA scatter + NJ tree) from a distance matrix --------------------
plot_level <- function(M, lvl, pca_path, nj_path) {
  n <- if (is.null(M)) 0 else nrow(M)
  if (n >= 3) {
    # PCoA (labels repelled + generous expansion so nothing runs off the page)
    co <- tryCatch(cmdscale(as.dist(M), k = 2, eig = TRUE), error = function(e) NULL)
    if (!is.null(co)) {
      pts <- as.data.frame(co$points); names(pts) <- c("Dim1", "Dim2"); pts$u <- short(rownames(pts))
      ev <- co$eig; ev[ev < 0] <- 0
      lab <- function(k) if (sum(ev) > 0) sprintf("PCoA %d (%.1f%%)", k, 100 * ev[k]/sum(ev)) else sprintf("PCoA %d", k)
      p <- ggplot(pts, aes(Dim1, Dim2, label = u)) +
        geom_point(size = 3, colour = "#2c7fb8") +
        geom_text_repel(size = 3.2, max.overlaps = Inf, min.segment.length = 0,
                        box.padding = 0.6, point.padding = 0.3, seed = 1) +
        scale_x_continuous(expand = expansion(mult = 0.22)) +
        scale_y_continuous(expand = expansion(mult = 0.14)) +
        labs(title = sprintf("%s %s %s ordination (PCoA)", label, DASH, lvl), x = lab(1), y = lab(2))
      ggsave(pca_path, p, width = 7.5, height = 5.5, dpi = 150)
    } else placeholder(pca_path, "PCoA failed")
    # NJ tree (short tips + margins + xpd so axial labels are not clipped)
    tr <- tryCatch(nj(as.dist(M)), error = function(e) NULL)
    if (!is.null(tr)) {
      tr$tip.label <- short(tr$tip.label)
      png(nj_path, width = 8.5 * 150, height = 6 * 150, res = 150)
      old <- par(mar = c(2, 2, 3, 2), xpd = NA)
      plot(tr, type = "unrooted", lab4ut = "axial", cex = 0.9, edge.width = 1.2, no.margin = FALSE,
           main = sprintf("%s %s %s NJ tree (graph similarity)", label, DASH, lvl), cex.main = 0.9)
      par(old); dev.off()
    } else placeholder(nj_path, "NJ tree failed")
  } else {
    placeholder(pca_path, sprintf("%s ordination needs %s3 units\n(have %d)", lvl, GE, n))
    placeholder(nj_path,  sprintf("%s NJ tree needs %s3 units\n(have %d)", lvl, GE, n))
  }
}

# ---- haplotype level -------------------------------------------------------------------
plot_level(D, "haplotype", op(".pca_haplotype.png"), op(".njtree_haplotype.png"))

# ---- individual level (aggregate haplotypes -> diploid individual) ---------------------
Di <- NULL
if (!is.null(D) && length(groups) >= 1) {
  individual_of <- function(h) sub("_(dip_)?hap[0-9]+$", "", sub("#.*$", "", h))
  inds <- individual_of(groups); uind <- sort(unique(inds))
  Di <- matrix(0, length(uind), length(uind), dimnames = list(uind, uind))
  for (i in seq_along(uind)) for (j in seq_along(uind)) if (i != j) {
    hi <- which(inds == uind[i]); hj <- which(inds == uind[j])
    Di[i, j] <- mean(D[hi, hj, drop = FALSE])
  }
}
plot_level(Di, "individual", op(".pca_individual.png"), op(".njtree_individual.png"))
