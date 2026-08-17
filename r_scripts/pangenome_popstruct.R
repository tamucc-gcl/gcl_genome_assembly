#!/usr/bin/env Rscript
# ======================================================================================
# pangenome_popstruct.R  (workstream D)
# Repo location: r_scripts/pangenome_popstruct.R
#
# From an `odgi similarity -D '#' -p 2 -d` table (grouped by PanSN haplotype), renders
# per-HAPLOTYPE and per-INDIVIDUAL ordination (PCoA) + midpoint-rooted neighbour-joining
# trees (rightwards phylogram, tip labels aligned). In the
# haplotype-level figures each unit is coloured by its diploid INDIVIDUAL and shaped by its
# HAPLOTYPE (both on the points/tips and in the labels); consistent colours are shared
# between the PCoA and the tree. Small-n tolerant (placeholder if < 3 units).
#
# Usage: Rscript pangenome_popstruct.R <similarity.tsv> <label> <outdir>
# ======================================================================================

suppressPackageStartupMessages({ library(ggplot2); library(ape); library(ggrepel) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: pangenome_popstruct.R <similarity.tsv> <label> <outdir>")
sim_f <- args[1]; label <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
op <- function(x) file.path(outdir, paste0(label, x))
theme_set(theme_bw(base_size = 12))
DASH <- "\u2014"; GE <- "\u2265"

# Short display labels, driven by the PanSN convention pangenome.nf writes into the cactus
# seqfile:
#   * a normal individual's haplotypes group as <individual>.<hap>  -> graph "<ind>#<N>";
#   * the REFERENCE individual's haplotypes cannot be grouped with the reference (cactus
#     forbids a sample sharing a name prefix with the reference), so they get FLAT names,
#     i.e. "<ind>_hap<N>" -> graph "<ind>_hap<N>#0". Collapse that back to "<ind>#<N>" so
#     the reference individual is not split into two pseudo-individuals carrying haplotype
#     0. The legacy "_dip_hap<N>#0" spelling is accepted too.
#   * an optional leading species tag ("Sde-") is dropped, but only when EVERY unit carries
#     the same tag -- no species is hard-coded.
strip_tag <- function(x) {
  if (!length(x)) return(x)
  tag <- sub("^([A-Za-z]{2,6}-).*$", "\\1", x)
  if (all(grepl("^[A-Za-z]{2,6}-", x)) && length(unique(tag)) == 1L)
    substring(x, nchar(tag[1]) + 1L) else x
}
short  <- function(x) sub("_(dip_)?hap([0-9]+)#0$", "#\\2", strip_tag(x))
ind_of <- function(u) sub("#.*$", "", u)                                  # individual
hap_of <- function(u) ifelse(grepl("#", u), sub("^.*#", "", u), NA_character_)  # "1"/"2"/NA
SHP    <- c("0" = 8, "1" = 16, "2" = 17, "3" = 15, "4" = 18)  # asterisk, circle, triangle, square, diamond
shp_of <- function(h) {                       # unmapped haplotype -> 'x'; no haplotype -> circle
  v <- unname(SHP[as.character(h)]); v[is.na(v)] <- 4L; v[is.na(h)] <- 16L; v
}

placeholder <- function(path, msg) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 5) +
    theme_void() + xlim(-1, 1) + ylim(-1, 1)
  ggsave(path, p, width = 7.8, height = 5.6, dpi = 150)
}

# ---- midpoint rooting (ape primitives only; no phangorn / phytools) --------------------
# Roots at the midpoint of the longest tip-to-tip path: the correct default for conspecific
# assemblies, where no outgroup exists. Total tree length is preserved. Returns the input
# unchanged if the tree cannot be rooted meaningfully (no/degenerate branch lengths, < 3
# tips, unexpected node path), so a failure degrades to the unrooted tree rather than
# killing the process.
midpoint_root <- function(tr) {
  if (is.null(tr$edge.length) || Ntip(tr) < 3L) return(tr)
  # collapse any existing degree-2 root first (nj() output is unrooted, but be safe): its
  # two edges are one edge unrooted, and treating them separately would corrupt the
  # root-edge length reassignment below. unroot() preserves total tree length.
  if (is.rooted(tr)) tr <- tryCatch(unroot(tr), error = function(e) tr)
  el <- tr$edge.length
  if (is.null(el)) return(tr)
  if (!all(is.finite(el)) || max(el) <= 0) return(tr)
  d <- tryCatch(as.matrix(cophenetic(tr)), error = function(e) NULL)
  if (is.null(d) || !all(is.finite(d)) || max(d) <= 0) return(tr)
  ij   <- which(d == max(d), arr.ind = TRUE)[1, ]
  half <- d[as.integer(ij[1]), as.integer(ij[2])] / 2
  i    <- as.integer(ij[1]); j <- as.integer(ij[2])
  if (!is.null(rownames(d))) {                 # index tips by label, not by row position
    mi <- match(rownames(d)[i], tr$tip.label); mj <- match(colnames(d)[j], tr$tip.label)
    if (is.na(mi) || is.na(mj)) return(tr)
    i <- mi; j <- mj
  }
  path <- tryCatch(as.integer(nodepath(tr, i, j)), error = function(e) NULL)
  if (is.null(path) || length(path) < 2L) return(tr)
  if (path[1] != i) { if (path[length(path)] == i) path <- rev(path) else return(tr) }
  cum <- 0
  for (k in seq_len(length(path) - 1L)) {
    e <- which((tr$edge[, 1] == path[k]      & tr$edge[, 2] == path[k + 1L]) |
               (tr$edge[, 1] == path[k + 1L] & tr$edge[, 2] == path[k]))
    if (length(e) != 1L) return(tr)
    L <- el[e]
    if (cum + L >= half) {
      lenA  <- half - cum                 # root .. path[k]   (the tip-i side of this edge)
      lenB  <- L - lenA                   # root .. path[k+1] (the tip-j side)
      child <- tr$edge[e, 2]              # ape's child end of e: the clade below it defines
      clen  <- if (child == path[k]) lenA else lenB      # one side of the new root exactly
      ctips <- if (child <= Ntip(tr)) tr$tip.label[child]
               else tryCatch(extract.clade(tr, child)$tip.label, error = function(x) NULL)
      if (is.null(ctips) || !length(ctips) || length(ctips) >= Ntip(tr)) return(tr)
      # root by outgroup = that clade: unambiguous, and (unlike node=) also legal when the
      # midpoint falls on a pendant edge, which is common for NJ trees
      rt <- tryCatch(suppressWarnings(root(tr, outgroup = ctips, resolve.root = TRUE)),
                     error = function(x) NULL)
      if (is.null(rt)) return(tr)
      rnode <- setdiff(unique(rt$edge[, 1]), rt$edge[, 2])
      re    <- which(rt$edge[, 1] == rnode[1])
      if (length(rnode) == 1L && length(re) == 2L) {
        n1 <- rt$edge[re[1], 2]
        t1 <- if (n1 <= Ntip(rt)) rt$tip.label[n1]
              else tryCatch(extract.clade(rt, n1)$tip.label, error = function(x) NULL)
        if (!is.null(t1)) {
          s <- if (setequal(t1, ctips)) c(1L, 2L) else c(2L, 1L)
          rt$edge.length[re[s[1]]] <- clen
          rt$edge.length[re[s[2]]] <- L - clen
        }
      }
      return(rt)
    }
    cum <- cum + L
  }
  tr
}

# ---- distance matrix from the odgi similarity table ------------------------------------
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

# ---- plot one level: PCoA (colour=individual, shape=haplotype) + NJ tree ---------------
plot_level <- function(M, lvl, pca_path, nj_path) {
  n <- if (is.null(M)) 0 else nrow(M)
  if (n < 3) {
    placeholder(pca_path, sprintf("%s ordination needs %s3 units\n(have %d)", lvl, GE, n))
    placeholder(nj_path,  sprintf("%s NJ tree needs %s3 units\n(have %d)", lvl, GE, n))
    return(invisible())
  }
  units <- short(rownames(M)); ind <- ind_of(units); hapn <- hap_of(units)
  ind_levels <- sort(unique(ind))
  col_pal <- setNames(scales::hue_pal()(length(ind_levels)), ind_levels)   # shared hue map
  has_hap <- any(!is.na(hapn))
  uh <- sort(unique(na.omit(hapn)))

  ## PCoA
  co <- tryCatch(cmdscale(as.dist(M), k = 2, eig = TRUE), error = function(e) NULL)
  if (!is.null(co)) {
    ev <- co$eig; ev[ev < 0] <- 0
    axl <- function(k) if (sum(ev) > 0) sprintf("PCoA %d (%.1f%%)", k, 100 * ev[k]/sum(ev)) else sprintf("PCoA %d", k)
    pts <- data.frame(Dim1 = co$points[, 1], Dim2 = co$points[, 2], u = units,
                      individual = factor(ind, ind_levels),
                      hap = ifelse(is.na(hapn), NA, paste0("hap", hapn)))
    p <- ggplot(pts, aes(Dim1, Dim2, colour = individual))
    if (has_hap) {
      p <- p + geom_point(aes(shape = hap), size = 3.8) +
        scale_shape_manual(values = setNames(shp_of(uh), paste0("hap", uh)), name = "haplotype")
    } else {
      p <- p + geom_point(size = 3.8, shape = 16)
    }
    p <- p +
      geom_text_repel(aes(label = u), size = 3.2, max.overlaps = Inf, min.segment.length = 0,
                      box.padding = 0.6, point.padding = 0.4, show.legend = FALSE, seed = 1) +
      scale_colour_manual(values = col_pal, name = "individual") +
      scale_x_continuous(expand = expansion(mult = 0.22)) +
      scale_y_continuous(expand = expansion(mult = 0.14)) +
      labs(title = sprintf("%s %s %s ordination (PCoA)", label, DASH, lvl), x = axl(1), y = axl(2))
    ggsave(pca_path, p, width = 7.8, height = 5.6, dpi = 150)
  } else placeholder(pca_path, "PCoA failed")

  ## NJ tree: midpoint-rooted rightwards phylogram, tip labels aligned on the right
  tr <- tryCatch(nj(as.dist(M)), error = function(e) NULL)
  if (!is.null(tr)) {
    tr$tip.label <- short(tr$tip.label)
    nneg <- sum(tr$edge.length < 0, na.rm = TRUE)     # nj can emit negative branch lengths;
    if (nneg > 0) {                                   # clamp so the rooting geometry is sane
      message(sprintf("[popstruct] %s NJ: clamped %d negative branch length(s) to 0", lvl, nneg))
      tr$edge.length[tr$edge.length < 0] <- 0
    }
    tr   <- tryCatch(ladderize(midpoint_root(tr), right = FALSE), error = function(e) tr)
    tind <- ind_of(tr$tip.label); thapn <- hap_of(tr$tip.label)
    tcol <- col_pal[tind]; tpch <- shp_of(thapn)
    dep  <- suppressWarnings(max(node.depth.edgelength(tr)))
    if (!is.finite(dep) || dep <= 0) dep <- 1
    png(nj_path, width = 8.5 * 150, height = 6 * 150, res = 150)
    old <- par(mar = c(5.5, 0.6, 3, 0.6), xpd = NA)
    plot(tr, type = "phylogram", direction = "rightwards", align.tip.label = TRUE,
         cex = 0.9, font = 1, tip.color = tcol, label.offset = 0.03 * dep, no.margin = FALSE,
         main = sprintf("%s %s %s NJ tree (graph similarity, midpoint-rooted)", label, DASH, lvl),
         cex.main = 0.9)
    tiplabels(pch = tpch, col = tcol, cex = 1.3)
    usr <- par("usr"); dy <- par("cxy")[2]
    add.scale.bar(x = usr[1], y = usr[4] + 0.9 * dy, length = signif(dep / 4, 1), cex = 0.8, lwd = 2)
    legend(x = usr[1], y = usr[3] - 1.6 * dy, legend = ind_levels, col = col_pal[ind_levels],
           pch = 16, pt.cex = 1.2, bty = "n", cex = 0.85, horiz = TRUE, title = "individual")
    if (has_hap) legend(x = usr[1], y = usr[3] - 3.4 * dy, legend = paste0("hap", uh),
                        pch = shp_of(uh), col = "grey30", pt.cex = 1.2, bty = "n", cex = 0.85,
                        horiz = TRUE, title = "haplotype")
    par(old); dev.off()
  } else placeholder(nj_path, "NJ tree failed")
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
