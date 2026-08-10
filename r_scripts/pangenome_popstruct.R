#!/usr/bin/env Rscript
# ======================================================================================
# pangenome_popstruct.R  (workstream D)
# Repo location: r_scripts/pangenome_popstruct.R
#
# From an `odgi similarity -D '#' -p 2 -d` long-format table (grouped by PanSN haplotype),
# renders per-HAPLOTYPE and per-INDIVIDUAL ordination (PCoA) + neighbour-joining trees.
#   - Haplotype level: each haploid assembly is one unit (incl. the reference).
#   - Individual level: haplotypes aggregated to their diploid individual (mean between-
#     individual pairwise distance). The individual is derived from the haplotype name by
#     stripping the PanSN hap markers (`#N`, and a reference-style `_dip_hapN` / `_hapN`).
# Distances come from shared graph content (SNPs + indels + SVs). Small-n tolerant: a level
# with < 3 units gets a labelled placeholder. Always writes all four PNGs.
#
# Usage: Rscript pangenome_popstruct.R <similarity.tsv> <label> <outdir>
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

# ---- build the haplotype x haplotype distance matrix -----------------------------------
D <- NULL; groups <- character(0)
if (file.exists(sim_f) && file.size(sim_f) > 0) {
  s <- tryCatch(read.delim(sim_f, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (!is.null(s) && nrow(s) > 0 && ncol(s) >= 3) {
    nm <- names(s)
    ga_c <- if ("group.a" %in% nm) "group.a" else nm[1]
    gb_c <- if ("group.b" %in% nm) "group.b" else nm[2]
    # distance column: prefer jaccard.distance (bounded 0-1), then euclidean.distance, then
    # any *.distance column. These are already DISTANCES (0 = identical), so use as-is.
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
    co <- tryCatch(cmdscale(as.dist(M), k = 2, eig = TRUE), error = function(e) NULL)
    if (!is.null(co)) {
      pts <- as.data.frame(co$points); names(pts) <- c("Dim1", "Dim2"); pts$u <- rownames(pts)
      ev <- co$eig; ev[ev < 0] <- 0
      lab <- function(k) if (sum(ev) > 0) sprintf("PCoA %d (%.1f%%)", k, 100 * ev[k]/sum(ev)) else sprintf("PCoA %d", k)
      p <- ggplot(pts, aes(Dim1, Dim2, label = u)) +
        geom_point(size = 3, colour = "#2c7fb8") + geom_text(vjust = -0.8, size = 3) +
        labs(title = sprintf("%s \u2014 %s ordination (PCoA)", label, lvl), x = lab(1), y = lab(2))
      ggsave(pca_path, p, width = 7, height = 5, dpi = 150)
    } else placeholder(pca_path, "PCoA failed")
    tr <- tryCatch(nj(as.dist(M)), error = function(e) NULL)
    if (!is.null(tr)) {
      png(nj_path, width = 7 * 150, height = 5 * 150, res = 150)
      plot(tr, type = "unrooted", lab4ut = "axial", main = sprintf("%s \u2014 %s NJ tree (graph similarity)", label, lvl))
      dev.off()
    } else placeholder(nj_path, "NJ tree failed")
  } else {
    placeholder(pca_path, sprintf("%s ordination needs \u22653 units\n(have %d)", lvl, n))
    placeholder(nj_path,  sprintf("%s NJ tree needs \u22653 units\n(have %d)", lvl, n))
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
