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
#   Rscript pangenome_plots.R <hist.tsv> <sv_sizes.tsv> <variant_summary.tsv> <label> <outdir> \
#       [hap_private=<tsv>] [ref=<sample#hap>] [core=1.0] [softcore=-1] [shell=2]
#
# Inputs:
#   hist.tsv            panacus coverage histogram: col1 = coverage level (1..N), col2 = count (bp)
#   sv_sizes.tsv        PANGENOME_VARIANTS: columns sv_size_bp, sv_type
#   variant_summary.tsv PANGENOME_VARIANTS: columns class, count
# Outputs (in outdir):
#   <label>.growth_curves.png, <label>.coverage_histogram.png,
#   <label>.sv_size_histogram.png, <label>.variant_summary.png, <label>.growth_fit.tsv
#   <label>.pangenome_partition.png, <label>.coverage_histogram_tiers.png
#   <label>.private_by_haplotype.png, <label>.private_fraction_by_haplotype.png
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

# ---- extra options: every arg after #5 is key=value (order-free, back-compatible) -----
#   hap_private=<tsv>       PANGENOME_HAP_COVERAGE per-haplotype private breakdown
#   ref=<sample#hap>        PanSN name of the reference haplotype (highlighted if given)
#   core= softcore= shell=  tier cuts (see cut_to_k below for the three accepted forms)
kv <- list()
for (a in args[-(1:5)]) if (grepl("=", a, fixed = TRUE)) kv[[sub("=.*$", "", a)]] <- sub("^[^=]*=", "", a)
optc <- function(k, d) if (!is.null(kv[[k]]) && nzchar(kv[[k]])) kv[[k]] else d
optn <- function(k, d) { v <- suppressWarnings(as.numeric(optc(k, NA))); if (length(v) == 1 && is.finite(v)) v else d }

hap_private_f <- optc("hap_private", "NONE")
ref_hap       <- optc("ref",         "NONE")
cut_core      <- optn("core",      1.00)
cut_softcore  <- optn("softcore",  -1)
cut_shell     <- optn("shell",      2)

# NO_* / NONE are how the pipeline signals "this sub-analysis did not run"
is_missing <- function(x) is.null(x) || is.na(x) || !nzchar(x) ||
  grepl("^(NONE|NO_)", basename(x)) || !file.exists(x) || file.size(x) == 0
mb <- function(x) paste0(formatC(x / 1e6, format = "f", digits = 0, big.mark = ","), " Mb")

tier_levels <- c("core", "soft-core", "shell", "cloud")
tier_cols   <- c(core = "#d95f0e", "soft-core" = "#fdae61", shell = "#c2a5cf", cloud = "#7fbf7b")

# A cut becomes a MINIMUM HAPLOTYPE COUNT k. Three accepted forms, so that the same knob
# works at n = 5 and at n = 500:
#   cut >  1   absolute count            (shell = 2  -> shared by at least two haplotypes)
#   0 < cut <= 1   fraction of N         (core = 1.0 -> all;  0.15 -> the classic cloud cut)
#   cut <  0   N - |cut|, "all but k"    (softcore = -1 -> all but one)
# The negative form exists because no fraction expresses "all but one" across cohort sizes:
# at N = 6, 0.9 rounds up to 6 and the soft-core tier silently vanishes into core.
cut_to_k <- function(cut, N) {
  k <- if (cut < 0) N + cut else if (cut <= 1) ceiling(cut * N - 1e-9) else ceiling(cut)
  as.integer(min(N, max(1, k)))
}
# clamped monotone so the ifelse cascade can never invert
tier_k <- function(N) {
  kc <- cut_to_k(cut_core, N)
  ks <- min(cut_to_k(cut_softcore, N), kc)
  kh <- min(cut_to_k(cut_shell, N), ks)
  c(core = kc, "soft-core" = ks, shell = kh)
}
tier_of <- function(i, N) {
  k <- tier_k(N)
  factor(ifelse(i >= k[["core"]],      "core",
         ifelse(i >= k[["soft-core"]], "soft-core",
         ifelse(i >= k[["shell"]],     "shell", "cloud"))), levels = tier_levels)
}

# ---- helpers: read a panacus TSV keying on data rows starting with an integer ----------
read_int_keyed <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[!grepl("^#", ln)]
  ln <- ln[grepl("^[0-9]+\t", ln)]                       # data rows: integer in col 1
  if (length(ln) == 0) return(NULL)
  do.call(rbind, lapply(strsplit(ln, "\t"), function(r) as.numeric(r[1:2])))
}

# ---- helper: '#'-commented TSV. comment.char is left OFF because PanSN haplotype names
# ---- contain '#' and would otherwise be truncated mid-field.
read_tsv_hash <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[!grepl("^#", ln)]
  if (length(ln) < 2) return(NULL)
  d <- try(read.delim(text = paste(ln, collapse = "\n"), header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE, comment.char = ""),
           silent = TRUE)
  if (inherits(d, "try-error") || nrow(d) == 0) NULL else d
}

# ======================================================================================
# 1. coverage histogram -> growth/core curves + Heaps fit + CI band
# ======================================================================================
h_mat <- read_int_keyed(hist_f)
if (!is.null(h_mat)) h_mat <- h_mat[h_mat[, 1] >= 1, , drop = FALSE]   # panacus emits a coverage-0 row; drop it (R is 1-indexed)
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

  # ---- core / soft-core / shell / cloud partition (settable cuts) ----------------------
  # Same h(i) as the plot above, re-binned by how many haplotypes carry a segment.
  # "pangenome" is the union of all four tiers, so it is not a fifth bin.
  tdf  <- data.frame(coverage = idx, bp = h, tier = tier_of(idx, N))
  sums <- tapply(tdf$bp, tdf$tier, sum)
  part <- data.frame(tier = factor(tier_levels, levels = tier_levels),
                     bp   = as.numeric(sums[tier_levels]))
  part$bp[is.na(part$bp)] <- 0
  part$pct <- 100 * part$bp / sum(part$bp)
  # floors are read back off the REALISED assignment, not recomputed from the cuts: at small
  # N a cut can be shadowed by a stricter one and the tier is then genuinely empty.
  realised <- function(tt) { lv <- idx[tdf$tier == tt]; if (length(lv)) min(lv) else NA_integer_ }
  floors   <- c(core = realised("core"), "soft-core" = realised("soft-core"),
                shell = realised("shell"))
  hapfmt  <- function(x) ifelse(is.na(x), "empty at this N", sprintf("\u2265%d/%d hap", x, N))
  thr_txt <- sprintf("core %s \u00b7 soft-core %s \u00b7 shell %s \u00b7 cloud below (cuts %g / %g / %g)",
                     hapfmt(floors[["core"]]), hapfmt(floors[["soft-core"]]),
                     hapfmt(floors[["shell"]]), cut_core, cut_softcore, cut_shell)

  p5 <- ggplot(part, aes(x = "", y = bp, fill = tier)) +
    geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(pct >= 4, sprintf("%s\n%.1f%%", mb(bp), pct), "")),
              position = position_stack(vjust = 0.5), size = 3.2, colour = "grey15") +
    coord_flip() +
    scale_y_continuous(labels = function(x) paste0(round(x / 1e6), " Mb"),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = tier_cols, drop = FALSE) +
    labs(title = paste0(label, " \u2014 pangenome partition"),
         subtitle = paste0(thr_txt, "\ntotal ", mb(sum(part$bp)), " over ", N, " haplotypes"),
         x = NULL, y = "sequence", fill = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())
  ggsave(op(".pangenome_partition.png"), p5, width = 8, height = 3.4, dpi = 150)

  # the same U-curve, coloured by tier, with the cuts drawn in so they are auditable
  vl <- as.numeric(floors[!is.na(floors)]) - 0.5
  p6 <- ggplot(tdf, aes(factor(coverage), bp, fill = tier)) +
    geom_col() +
    (if (length(vl)) geom_vline(xintercept = vl, linetype = 2, colour = "grey35") else NULL) +
    scale_y_continuous(labels = function(x) paste0(round(x / 1e6), " Mb")) +
    scale_fill_manual(values = tier_cols, drop = FALSE) +
    labs(title = paste0(label, " \u2014 coverage histogram by tier"), subtitle = thr_txt,
         x = "number of haplotypes covering a segment", y = "sequence", fill = NULL)
  ggsave(op(".coverage_histogram_tiers.png"), p6, width = 7, height = 5, dpi = 150)

  stopifnot(length(part$bp) == length(tier_levels))
  growth_fit <- rbind(growth_fit, data.frame(
    metric = c("tier_core_bp", "tier_softcore_bp", "tier_shell_bp", "tier_cloud_bp",
               "tier_core_cut", "tier_softcore_cut", "tier_shell_cut",
               "tier_core_min_haps", "tier_softcore_min_haps", "tier_shell_min_haps"),
    value  = c(part$bp, cut_core, cut_softcore, cut_shell,
               floors[["core"]], floors[["soft-core"]], floors[["shell"]])))
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

# ======================================================================================
# 4. private-sequence ownership by haplotype  (from PANGENOME_HAP_COVERAGE)
# ======================================================================================
# panacus's h(1) says HOW MUCH sequence is private but not WHOSE it is. gfa_hap_coverage.py
# resolves the owner per node, so these two plots decompose that single bar.
if (!is_missing(hap_private_f)) {
  hp   <- read_tsv_hash(hap_private_f)
  need <- c("haplotype", "private_bp", "hap_bp", "pct_of_private", "pct_of_haplotype")
  if (!is.null(hp) && all(need %in% names(hp)) && nrow(hp) > 0) {
    for (cc in need[-1]) hp[[cc]] <- suppressWarnings(as.numeric(hp[[cc]]))
    hp$name   <- as.character(hp$haplotype)
    hp$is_ref <- hp$name == ref_hap
    nh        <- nrow(hp)
    hgt       <- max(3.4, 0.42 * nh + 1.5)
    fill_sc   <- scale_fill_manual(values = c("FALSE" = "#7fbf7b", "TRUE" = "#2c7fb8"),
                                   labels = c("FALSE" = "haplotype", "TRUE" = "reference"),
                                   guide  = if (any(hp$is_ref)) "legend" else "none")

    # (a) share of the pangenome's private sequence. Every private segment has exactly one
    #     owner, so these bars sum to 100% by construction.
    hp$y1 <- factor(hp$name, levels = hp$name[order(hp$pct_of_private)])
    even  <- 100 / nh
    p7 <- ggplot(hp, aes(pct_of_private, y1, fill = is_ref)) +
      geom_col(width = 0.72) + fill_sc +
      geom_vline(xintercept = even, linetype = 2, colour = "grey35") +
      geom_text(aes(label = mb(private_bp)), hjust = -0.12, size = 3.1, colour = "grey20") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
      labs(title = paste0(label, " \u2014 who owns the private sequence"),
           subtitle = sprintf("total private %s; dashed line = even share (%.1f%%)",
                              mb(sum(hp$private_bp, na.rm = TRUE)), even),
           x = "% of the pangenome's private sequence", y = NULL, fill = NULL)
    ggsave(op(".private_by_haplotype.png"), p7, width = 8, height = hgt, dpi = 150)

    # (b) the same bp normalised BY HAPLOTYPE: how much of each assembly's graph content is
    #     unique to it. Unlike (a) this does not depend on the other haplotypes' sizes.
    hp$y2 <- factor(hp$name, levels = hp$name[order(hp$pct_of_haplotype)])
    p8 <- ggplot(hp, aes(pct_of_haplotype, y2, fill = is_ref)) +
      geom_col(width = 0.72) + fill_sc +
      geom_text(aes(label = sprintf("%s / %s", mb(private_bp), mb(hap_bp))),
                hjust = -0.08, size = 3.0, colour = "grey20") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.26))) +
      labs(title = paste0(label, " \u2014 private fraction of each haplotype"),
           subtitle = "private bp / this haplotype's total graph bp",
           x = "% of the haplotype that is private to it", y = NULL, fill = NULL)
    ggsave(op(".private_fraction_by_haplotype.png"), p8, width = 8, height = hgt, dpi = 150)
  } else message("hap_private table unusable; skipping private-by-haplotype plots")
} else message("no hap_private table; skipping private-by-haplotype plots")
