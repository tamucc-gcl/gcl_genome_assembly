#!/usr/bin/env Rscript
# ======================================================================================
# pangenome_plots.R  (workstream D)
# Repo location: r_scripts/pangenome_plots.R
#
# Turns the panacus + variant-catalog tables into the pangenome report figures. The growth
# and core curves, the Heaps'-law fit, and the confidence band are all derived from the
# coverage histogram h(i) = sequence covered by exactly i haplotypes, via the standard
# rarefaction formulas (this matches panacus's exact expected growth and adds a Coleman
# variance band, which panacus does not emit). Robust to small n.
#
# Usage:
#   Rscript pangenome_plots.R <hist.tsv> <sv_sizes.tsv> <variant_summary.tsv> <label> <outdir>
#
# Inputs:
#   hist.tsv            panacus coverage histogram: col1 = coverage level (1..N), col2 = count (bp)
#   sv_sizes.tsv        PANGENOME_VARIANTS: columns sv_size_bp, sv_type
#   variant_summary.tsv PANGENOME_VARIANTS: columns class, count
# Outputs (in outdir):
#   <label>.growth_curves.pdf, <label>.coverage_histogram.pdf,
#   <label>.sv_size_histogram.pdf, <label>.variant_summary.pdf, <label>.growth_fit.tsv
# ======================================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

args   <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: pangenome_plots.R <hist.tsv> <sv_sizes.tsv> <variant_summary.tsv> <label> <outdir>")
hist_f <- args[1]; sv_f <- args[2]; vs_f <- args[3]; label <- args[4]; outdir <- args[5]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
op <- function(x) file.path(outdir, paste0(label, x))

theme_set(theme_bw(base_size = 12))

# ---- helpers: read a panacus TSV keying on data rows starting with an integer ----------
read_int_keyed <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[!grepl("^#", ln)]
  ln <- ln[grepl("^[0-9]+\t", ln)]                       # data rows: integer in col 1
  if (length(ln) == 0) return(NULL)
  do.call(rbind, lapply(strsplit(ln, "\t"), function(r) as.numeric(r[1:2])))
}

# ======================================================================================
# 1. coverage histogram -> growth/core curves + Heaps fit + CI band
# ======================================================================================
h_mat <- read_int_keyed(hist_f)
growth_fit <- data.frame()

if (!is.null(h_mat) && nrow(h_mat) >= 1) {
  N   <- max(h_mat[, 1])
  h   <- numeric(N); h[h_mat[, 1]] <- h_mat[, 2]         # h[i] = bp covered by exactly i haplotypes
  idx <- 1:N

  # rarefaction expectations over a random m-subset of N haplotypes (exact), + Coleman variance
  lch <- function(a, b) ifelse(b < 0 | b > a, -Inf, lchoose(a, b))
  pan_mean <- core_mean <- pan_sd <- numeric(N)
  for (m in 1:N) {
    p_cov <- 1 - exp(lch(N - idx, m) - lch(N, m))        # P(bp covered by >=1 of m)  per coverage class i
    p_all <-     exp(lch(idx,     m) - lch(N, m))        # P(bp in ALL m)             per coverage class i
    pan_mean[m]  <- sum(h * p_cov)
    core_mean[m] <- sum(h * p_all)
    pan_sd[m]    <- sqrt(sum(h * p_cov * (1 - p_cov)))    # Coleman (independence) band
  }

  # Heaps' law fit: pangenome size ~ k * m^gamma  (log-log least squares; needs m>=2 points)
  gamma <- NA_real_; kcoef <- NA_real_
  if (N >= 3) {
    fit <- try(lm(log(pan_mean) ~ log(idx)), silent = TRUE)
    if (!inherits(fit, "try-error")) { gamma <- unname(coef(fit)[2]); kcoef <- exp(unname(coef(fit)[1])) }
  }

  gdf <- rbind(
    data.frame(m = idx, value = pan_mean,  lo = pan_mean - 1.96 * pan_sd,
               hi = pan_mean + 1.96 * pan_sd, curve = "pangenome"),
    data.frame(m = idx, value = core_mean, lo = core_mean, hi = core_mean, curve = "core")
  )
  sub <- if (!is.na(gamma)) sprintf("Heaps: pan(m) \u2248 %.3g \u00b7 m^%.3f", kcoef, gamma) else "Heaps fit needs \u22653 haplotypes"
  p1 <- ggplot(gdf, aes(m, value, color = curve, fill = curve)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
    scale_x_continuous(breaks = idx) +
    scale_y_continuous(labels = function(x) paste0(round(x / 1e6), " Mb")) +
    scale_color_manual(values = c(pangenome = "#2c7fb8", core = "#d95f0e")) +
    scale_fill_manual(values  = c(pangenome = "#2c7fb8", core = "#d95f0e")) +
    labs(title = paste0(label, " \u2014 pangenome growth & core"),
         subtitle = sub, x = "haplotypes (m)", y = "sequence", color = NULL, fill = NULL)
  ggsave(op(".growth_curves.png"), p1, width = 7, height = 5, dpi = 150)

  # coverage histogram (U-curve): bp per haplotype-coverage level
  private <- h[1]; core <- h[N]; accessory <- if (N > 2) sum(h[2:(N - 1)]) else 0
  cls <- ifelse(idx == 1, "private", ifelse(idx == N, "core", "accessory"))
  hdf <- data.frame(coverage = idx, bp = h, class = cls)
  p2 <- ggplot(hdf, aes(factor(coverage), bp, fill = class)) +
    geom_col() +
    scale_y_continuous(labels = function(x) paste0(round(x / 1e6), " Mb")) +
    scale_fill_manual(values = c(private = "#7fbf7b", accessory = "#c2a5cf", core = "#d95f0e")) +
    labs(title = paste0(label, " \u2014 coverage histogram"),
         x = "number of haplotypes covering a segment", y = "sequence", fill = NULL)
  ggsave(op(".coverage_histogram.png"), p2, width = 7, height = 5, dpi = 150)

  growth_fit <- data.frame(
    metric = c("n_haplotypes", "pangenome_bp", "core_bp", "accessory_bp", "private_bp",
               "heaps_gamma", "heaps_k", "openness"),
    value  = c(N, pan_mean[N], core, accessory, private, gamma, kcoef,
               ifelse(is.na(gamma), "NA", ifelse(gamma > 0.1, "open", "closed"))))
} else {
  message("hist table empty/unreadable; skipping growth + coverage plots")
}
write.table(growth_fit, op(".growth_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ======================================================================================
# 2. SV size histogram (from sv_sizes.tsv: sv_size_bp, sv_type)
# ======================================================================================
sv <- try(read.delim(sv_f, stringsAsFactors = FALSE), silent = TRUE)
if (!inherits(sv, "try-error") && nrow(sv) > 0 && "sv_size_bp" %in% names(sv)) {
  sv$sv_size_bp <- suppressWarnings(as.numeric(sv$sv_size_bp))
  sv <- sv[is.finite(sv$sv_size_bp) & sv$sv_size_bp > 0, , drop = FALSE]
  if (nrow(sv) > 0) {
    p3 <- ggplot(sv, aes(sv_size_bp, fill = sv_type)) +
      geom_histogram(bins = 60, position = "stack", colour = NA) +
      scale_x_log10(labels = scales::comma) +
      labs(title = paste0(label, " \u2014 SV size spectrum"),
           x = "SV size (bp, log scale)", y = "count", fill = NULL)
    ggsave(op(".sv_size_histogram.png"), p3, width = 7, height = 5, dpi = 150)
  }
} else message("sv_sizes table empty/unreadable; skipping SV histogram")

# ======================================================================================
# 3. variant class summary bar (from variant_summary.tsv: class, count)
# ======================================================================================
vs <- try(read.delim(vs_f, stringsAsFactors = FALSE), silent = TRUE)
if (!inherits(vs, "try-error") && nrow(vs) > 0 && all(c("class", "count") %in% names(vs))) {
  keep <- c("SNP", "INDEL", "SV")
  vv <- vs[vs$class %in% keep, , drop = FALSE]
  vv$class <- factor(vv$class, levels = keep)
  vv$count <- as.numeric(vv$count)
  p4 <- ggplot(vv, aes(class, count, fill = class)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = scales::comma(count)), vjust = -0.3, size = 3.5) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
    labs(title = paste0(label, " \u2014 variant catalog"), x = NULL, y = "count")
  ggsave(op(".variant_summary.png"), p4, width = 6, height = 5, dpi = 150)
} else message("variant_summary table empty/unreadable; skipping variant bar")
