#!/usr/bin/env Rscript
# ======================================================================================
# pangenome_private_plots.R
# Repo location: r_scripts/pangenome_private_plots.R
#
# The private-sequence figure set. Takes BOTH graph flavours in one invocation, because the
# single most important thing these figures have to convey is that the two arms disagree:
# clipping removes 704,499,041 bp of which 698,360,436 -- 99.1% -- is private sequence, so the
# clip arm understates private content by 46% AND inverts the reference's apparent rank
# (highest of ten haplotypes on clip at 15.08%, lowest of ten on full at 8.08%, because the
# reference is the graph backbone and is never clipped). A figure drawn from one arm alone
# cannot show that.
#
# Usage:
#   Rscript pangenome_private_plots.R <label> <outdir> \
#       spectrum_clip=<tsv> spectrum_full=<tsv> \
#       hap_clip=<tsv>      hap_full=<tsv> \
#       contig_clip=<tsv>   contig_full=<tsv> \
#       evidence_clip=<csv> evidence_full=<csv> \
#       xtab_clip=<tsv>     xtab_full=<tsv> \
#       [ref=<sample#hap>] [bins=<csv of edges>] [min_bp=1000]
#
# Every input is optional and signalled missing the same way as pangenome_plots.R -- a NONE /
# NO_* basename, an absent file, or a zero-byte file. A figure whose inputs are missing is
# skipped with a message rather than failing the task.
#
# INPUT SCHEMAS (verified against real output, not assumed)
#   private_segment_spectrum.tsv  scope  size_bin  n_segments  segment_bp
#                                 scope is ALL or a haplotype key
#   hap_private.tsv               haplotype  sample  hap  private_bp  hap_bp
#                                 pct_of_private  pct_of_haplotype
#   hap_private_by_contig.tsv     haplotype  contig  private_bp  contig_graph_bp
#                                 pct_of_contig  pct_of_hap_private  pct_of_pangenome_private
#   private_evidence.csv          27 columns; see PANGENOME_PRIVATE_JOIN
#   private_evidence_xtab.tsv     scope  key  set  combined  n_segments  segment_bp
#
# Outputs (in outdir):
#   <label>.private_spectrum_bp.png            <label>.private_spectrum_count.png
#   <label>.private_by_chromosome.png          <label>.private_clip_vs_full.png
#   <label>.private_vs_control.png             <label>.private_evidence_xtab.png
#   <label>.private_by_haplotype.png           <label>.private_fraction_by_haplotype.png
#   <label>.private_figures_audit.tsv
#
# The last two MOVED here from pangenome_plots.R, which is pinned to the clip arm by its join
# with PANGENOME_GROWTH -- so they were being drawn from the arm that understates private
# sequence by 46% and inverts the reference's rank, which is what one of them is about.
# ======================================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: pangenome_private_plots.R <label> <outdir> key=value ...")
label <- args[1]; outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
op <- function(x) file.path(outdir, paste0(label, x))
theme_set(theme_bw(base_size = 12))

kv <- list()
for (a in args[-(1:2)]) if (grepl("=", a, fixed = TRUE))
  kv[[sub("=.*$", "", a)]] <- sub("^[^=]*=", "", a)
optc <- function(k, d = "NONE") if (!is.null(kv[[k]]) && nzchar(kv[[k]])) kv[[k]] else d
optn <- function(k, d) { v <- suppressWarnings(as.numeric(optc(k, NA))); if (length(v) == 1 && is.finite(v)) v else d }

# NO_* / NONE are how the pipeline signals "this sub-analysis did not run"
is_missing <- function(x) is.null(x) || is.na(x) || !nzchar(x) ||
  grepl("^(NONE|NO_)", basename(x)) || !file.exists(x) || file.size(x) == 0

# read.delim does NOT guarantee numeric types: with check.names = FALSE, one stray
# non-numeric token turns a whole column character, and arithmetic then fails at plot time.
# PANGENOME_REARRANGE_PLOTS died exactly this way ("'x' must be numeric" from cut()), so
# coerce explicitly on read rather than trusting the type.
num_cols <- function(d, cols) {
  for (cc in intersect(cols, names(d))) d[[cc]] <- suppressWarnings(as.numeric(d[[cc]]))
  d
}

rd <- function(f, sep = "\t", numeric_cols = character(0)) {
  if (is_missing(f)) return(NULL)
  d <- try(read.delim(f, sep = sep, header = TRUE, comment.char = "#",
                      stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
  if (inherits(d, "try-error") || is.null(d) || !nrow(d)) return(NULL)
  if (length(numeric_cols)) d <- num_cols(d, numeric_cols)
  d
}

ref_hap <- optc("ref", "NONE")
min_bp  <- optn("min_bp", 1000)
mb  <- function(x) paste0(formatC(x / 1e6, format = "f", digits = 0, big.mark = ","), " Mb")
gbp <- function(x) ifelse(x >= 1e9,
                          paste0(formatC(x / 1e9, format = "f", digits = 2), " Gb"),
                          paste0(formatC(x / 1e6, format = "f", digits = 0), " Mb"))

audit <- list()
note  <- function(k, v) audit[[k]] <<- as.character(v)
skip  <- function(fig, why) { message(sprintf("[private_plots] SKIP %s -- %s", fig, why))
                              note(paste0("skipped_", fig), why) }

# size_bin arrives as a string like "1000-5000" or ">=1000000"; order it numerically by its
# lower edge rather than alphabetically, or 100-500 sorts before 0-100.
bin_order <- function(x) {
  lo <- suppressWarnings(as.numeric(sub("^>=", "", sub("-.*$", "", x))))
  factor(x, levels = unique(x[order(lo)]))
}

# ---- both flavours, long ------------------------------------------------------------
bind_flavours <- function(clip_f, full_f, reader = rd) {
  parts <- list()
  for (fl in c("clip", "full")) {
    d <- reader(if (fl == "clip") clip_f else full_f)
    if (!is.null(d)) { d$flavor <- fl; parts[[fl]] <- d }
  }
  if (!length(parts)) return(NULL)
  do.call(rbind, parts)
}

spec <- bind_flavours(optc("spectrum_clip"), optc("spectrum_full"),
        reader = function(f) rd(f, numeric_cols = c("n_segments", "segment_bp")))
hap  <- bind_flavours(optc("hap_clip"), optc("hap_full"),
        reader = function(f) rd(f, numeric_cols = c("private_bp", "hap_bp",
                                "pct_of_private", "pct_of_haplotype")))
ctg  <- bind_flavours(optc("contig_clip"), optc("contig_full"),
        reader = function(f) rd(f, numeric_cols = c("private_bp", "contig_graph_bp",
                                "pct_of_contig", "pct_of_hap_private",
                                "pct_of_pangenome_private")))
xtab <- bind_flavours(optc("xtab_clip"), optc("xtab_full"),
        reader = function(f) rd(f, numeric_cols = c("n_segments", "segment_bp")))
ev   <- bind_flavours(optc("evidence_clip"), optc("evidence_full"),
        reader = function(f) rd(f, sep = ",",
                               numeric_cols = c("span_bp", "is_private", "repeat_like",
                                                "median_copy", "copy_ratio",
                                                "aligned_frac_merged",
                                                "n_other_assemblies")))

note("flavours_spectrum", if (is.null(spec)) 0 else length(unique(spec$flavor)))
note("flavours_evidence", if (is.null(ev))   0 else length(unique(ev$flavor)))

# ======================================================================================
# 1-2. The private-segment size spectrum, bp-weighted AND by count
#
# THE bp-WEIGHTED ONE IS THE DELIVERABLE. 99.0% of private segments are under 100 bp and
# carry 7% of the sequence, while the 1-50 kb range holds 601 Mb -- 74% of all private
# sequence -- in 0.5% of segments. A count histogram is a single spike at the left edge and
# says nothing. Both are drawn so that point is explicit rather than implied.
#
# Bin edges are pangenome_sv_bins, SHARED with the SV size spectrum, which is what makes
# "is the private-haplotype spectrum shaped like the SV spectrum?" readable off one axis.
# ======================================================================================
if (is.null(spec)) {
  skip("private_spectrum", "no private_segment_spectrum.tsv for either flavour")
} else {
  s <- spec[spec$scope == "ALL", ]
  if (!nrow(s)) {
    skip("private_spectrum", "no ALL-scope rows")
  } else {
    s <- s[is.finite(s$segment_bp) & is.finite(s$n_segments), , drop = FALSE]
    # scale_y_log10 silently drops n_segments <= 0; do it explicitly and record the count
    note("spectrum_bins_zero_count", sum(s$n_segments <= 0))
    s$size_bin <- bin_order(s$size_bin)
    note("spectrum_total_bp_clip", sum(s$segment_bp[s$flavor == "clip"]))
    note("spectrum_total_bp_full", sum(s$segment_bp[s$flavor == "full"]))

    p1 <- ggplot(s, aes(size_bin, segment_bp, fill = flavor)) +
      geom_col(position = position_dodge(preserve = "single")) +
      scale_y_continuous(labels = function(v) gbp(v)) +
      scale_fill_manual(values = c(clip = "grey55", full = "#2c7fb8")) +
      labs(title = paste0(label, " — private sequence by segment size (bp-weighted)"),
           subtitle = paste0("clipping removes 99.1% private sequence, so the clip arm ",
                             "understates by ~46%; floor ", min_bp, " bp for downstream sets"),
           x = "segment size (bp)", y = "private sequence", fill = "graph") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top")
    ggsave(op(".private_spectrum_bp.png"), p1, width = 9, height = 5, dpi = 150)

    p2 <- ggplot(s, aes(size_bin, n_segments, fill = flavor)) +
      geom_col(position = position_dodge(preserve = "single")) +
      scale_y_log10(labels = function(v) formatC(v, format = "d", big.mark = ",")) +
      scale_fill_manual(values = c(clip = "grey55", full = "#2c7fb8")) +
      labs(title = paste0(label, " — private segments by size (COUNT, log scale)"),
           subtitle = paste0("shown for contrast: 99.0% of segments are <100 bp and carry ",
                             "7% of the sequence. Counts answer a different question."),
           x = "segment size (bp)", y = "segments (log)", fill = "graph") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top")
    ggsave(op(".private_spectrum_count.png"), p2, width = 9, height = 5, dpi = 150)
  }
}

# ======================================================================================
# 3. Private fraction per chromosome
#
# The live question: chr8 and chr9 carry an elevated private fraction against a floor across
# the other chromosomes, consistently across all ten haplotypes and both flavours. Plotting
# pct_of_contig (private bp / contig graph bp) rather than raw bp is what makes it a fraction
# question rather than a chromosome-length question.
# ======================================================================================
if (is.null(ctg)) {
  skip("private_by_chromosome", "no hap_private_by_contig.tsv for either flavour")
} else {
  c2 <- ctg
  # placed chromosomes only: unplaced scaffolds are individually tiny and would swamp the
  # axis with hundreds of categories
  c2 <- c2[grepl("^chr", c2$contig) & is.finite(c2$pct_of_contig), ]
  if (!nrow(c2)) {
    skip("private_by_chromosome", "no chr* contigs")
  } else {
    # chr10_1 / chr10_17+chr11_12 -> chr10, so pieces of one chromosome pool
    c2$chrom <- sub("_.*$", "", c2$contig)
    ord <- unique(c2$chrom)
    num <- suppressWarnings(as.numeric(sub("^chr", "", ord)))
    c2$chrom <- factor(c2$chrom, levels = ord[order(num, ord)])

    p3 <- ggplot(c2, aes(chrom, pct_of_contig, colour = flavor)) +
      geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.75)) +
      scale_colour_manual(values = c(clip = "grey40", full = "#2c7fb8")) +
      labs(title = paste0(label, " — private fraction per chromosome"),
           subtitle = "one point per haplotype x chromosome piece; private bp / contig graph bp",
           x = NULL, y = "private fraction of contig (%)", colour = "graph") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
    ggsave(op(".private_by_chromosome.png"), p3, width = 9, height = 5, dpi = 150)

    agg <- aggregate(pct_of_contig ~ chrom + flavor, data = c2, FUN = median)
    note("chrom_median_spread_full",
         sprintf("%.2f-%.2f", min(agg$pct_of_contig[agg$flavor == "full"]),
                 max(agg$pct_of_contig[agg$flavor == "full"])))
  }
}

# ======================================================================================
# 4. Clip vs full, per haplotype -- the figure that shows the arms disagree
#
# Two panels: absolute private bp, and share of the cohort's private sequence. The second is
# where the reference rank inversion appears, and it is the reason the reference is
# highlighted when `ref=` is supplied.
# ======================================================================================
if (is.null(hap) || length(unique(hap$flavor)) < 2) {
  skip("private_clip_vs_full",
       if (is.null(hap)) "no hap_private.tsv" else "only one flavour available")
} else {
  h <- hap
  h$is_ref <- if (ref_hap == "NONE") FALSE else h$haplotype == ref_hap
  # order by full-arm private bp so the panels share a haplotype order
  h <- h[is.finite(h$private_bp), , drop = FALSE]
  o <- if ("full" %in% h$flavor) h[h$flavor == "full", ] else h
  h$haplotype <- factor(h$haplotype, levels = o$haplotype[order(o$private_bp)])

  p4 <- ggplot(h, aes(private_bp, haplotype, fill = flavor)) +
    geom_col(position = position_dodge(preserve = "single")) +
    geom_point(data = h[h$is_ref, ], aes(x = private_bp, y = haplotype),
               shape = 8, size = 2, colour = "firebrick", inherit.aes = FALSE) +
    scale_x_continuous(labels = function(v) gbp(v)) +
    scale_fill_manual(values = c(clip = "grey55", full = "#2c7fb8")) +
    labs(title = paste0(label, " — private sequence per haplotype, clip vs full"),
         subtitle = paste0("star = reference. Clipping never touches the reference path, ",
                           "which is why its RANK differs between arms."),
         x = "private sequence", y = NULL, fill = "graph") +
    theme(legend.position = "top")
  ggsave(op(".private_clip_vs_full.png"), p4, width = 9,
         height = max(4, 0.42 * length(unique(h$haplotype)) + 2.2), dpi = 150)

  for (fl in c("clip", "full")) {
    hh <- h[h$flavor == fl, ]
    if (nrow(hh) && ref_hap != "NONE" && any(hh$haplotype == ref_hap)) {
      r <- rank(-hh$pct_of_private)[hh$haplotype == ref_hap]
      note(paste0("reference_rank_", fl), sprintf("%d of %d", r, nrow(hh)))
    }
  }
}

# ======================================================================================
# 5. Private vs control -- the contrast the control set exists for
#
# Two independent measures side by side: k-mer repeat_like (copy number in the sample's own
# reads, relative to the control's single-copy level) and map NOT_PRIVATE (does the sequence
# align to another assembly). The control is size-matched, cross-individual, non-private
# sequence from the SAME haplotype, so any difference is attributable to privateness.
#
# NOTE the control's mapping ceiling. Cross-individual sharing in the graph does not imply
# block alignability: a window can be 95% covered by cross-individual nodes while every node
# is a few hundred bp, so it has no contiguous homologue. Read the private figure against the
# control's own value, not against 0 or 1.
# ======================================================================================
if (is.null(ev)) {
  skip("private_vs_control", "no private_evidence.csv for either flavour")
} else {
  e <- ev
  e$repeat_like_n <- suppressWarnings(as.numeric(e$repeat_like))
  rows <- list()
  for (fl in unique(e$flavor)) for (st in unique(e$set)) {
    sub <- e[e$flavor == fl & e$set == st, ]
    if (!nrow(sub)) next
    rl <- sub$repeat_like_n[is.finite(sub$repeat_like_n)]
    rows[[length(rows) + 1]] <- data.frame(
      flavor = fl, set = st, n = nrow(sub),
      measure = "k-mer REPEAT_LIKE",
      value = if (length(rl)) mean(rl) else NA_real_, stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      flavor = fl, set = st, n = nrow(sub),
      measure = "map NOT_PRIVATE",
      value = mean(sub$map_verdict == "NOT_PRIVATE"), stringsAsFactors = FALSE)
  }
  cmp <- do.call(rbind, rows)
  if (is.null(cmp) || !nrow(cmp)) {
    skip("private_vs_control", "no private/control rows")
  } else {
    p5 <- ggplot(cmp, aes(set, value, fill = set)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = sprintf("%.3f\n(n=%s)", value, formatC(n, big.mark = ","))),
                vjust = -0.15, size = 3) +
      facet_grid(measure ~ flavor) +
      scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
      scale_fill_manual(values = c(private = "#d95f0e", control = "grey60")) +
      labs(title = paste0(label, " — private vs matched non-private control"),
           subtitle = paste0("control = size-matched, cross-individual, non-private windows ",
                             "from the SAME haplotype, measured identically"),
           x = NULL, y = "fraction of segments") +
      theme(legend.position = "none")
    ggsave(op(".private_vs_control.png"), p5, width = 8, height = 6, dpi = 150)

    for (i in seq_len(nrow(cmp)))
      note(sprintf("%s_%s_%s", gsub("[^a-z]", "", tolower(cmp$measure[i])),
                   cmp$flavor[i], cmp$set[i]), sprintf("%.4f", cmp$value[i]))
  }
}

# ======================================================================================
# 6. The map x kmer verdict cross-tabulation, per chromosome
#
# The four combinations each mean something specific and only the cross-tab shows them:
#   NOT_PRIVATE       + REPEAT_LIKE   present elsewhere AND high copy -> graph collapse
#   PRIVATE_CONFIRMED + UNIQUE_LIKE   absent elsewhere AND single copy -> novel sequence
#   NOT_PRIVATE       + UNIQUE_LIKE   the graph failed to merge homologous sequence
#   PRIVATE_CONFIRMED + REPEAT_LIKE   haplotype-specific expansion
# The off-diagonals are the informative cases and are invisible in either measure alone.
# ======================================================================================
if (is.null(xtab)) {
  skip("private_evidence_xtab", "no private_evidence_xtab.tsv for either flavour")
} else {
  x <- xtab[xtab$scope == "CHROM" & xtab$set == "private" &
            is.finite(xtab$n_segments), ]
  if (!nrow(x)) {
    skip("private_evidence_xtab", "no CHROM-scope private rows")
  } else {
    tot <- aggregate(n_segments ~ key + flavor, data = x, FUN = sum)
    names(tot)[names(tot) == "n_segments"] <- "tot"
    x <- merge(x, tot, by = c("key", "flavor"))
    x$frac <- x$n_segments / x$tot
    ord <- unique(x$key)
    num <- suppressWarnings(as.numeric(sub("^chr", "", ord)))
    x$key <- factor(x$key, levels = ord[order(num, ord)])

    p6 <- ggplot(x, aes(key, frac, fill = combined)) +
      geom_col() +
      facet_wrap(~ flavor, ncol = 1) +
      scale_y_continuous(labels = function(v) sprintf("%.0f%%", 100 * v)) +
      labs(title = paste0(label, " — private-segment evidence by chromosome"),
           subtitle = "map verdict x k-mer verdict; off-diagonals are the informative cases",
           x = NULL, y = "share of private segments", fill = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right", legend.text = element_text(size = 8))
    ggsave(op(".private_evidence_xtab.png"), p6, width = 10, height = 7, dpi = 150)
  }
}

# ======================================================================================
# 7-8. Private ownership and private fraction per haplotype
#
# MOVED HERE FROM pangenome_plots.R, which is pinned to the clip arm because it is joined with
# PANGENOME_GROWTH (panacus on the clip GFA). These two figures were therefore being drawn from
# the arm that understates private sequence by 46% and inverts the reference's rank -- which is
# precisely what figure (a) is about. Faceted by flavour so the inversion is visible.
#
# (a) share of the cohort's private sequence. Every private segment has exactly one owner, so
#     the bars sum to 100% by construction, and the even-share line is the null.
# (b) the same bp normalised BY HAPLOTYPE: how much of each assembly's graph content is unique
#     to it. Unlike (a) this does not depend on the other haplotypes' sizes.
# ======================================================================================
if (is.null(hap)) {
  skip("private_ownership", "no hap_private.tsv for either flavour")
} else {
  need <- c("haplotype", "private_bp", "hap_bp", "pct_of_private", "pct_of_haplotype")
  if (!all(need %in% names(hap))) {
    skip("private_ownership", paste("missing columns:",
         paste(setdiff(need, names(hap)), collapse = ",")))
  } else {
    hp <- hap
    for (cc in need[-1]) hp[[cc]] <- suppressWarnings(as.numeric(hp[[cc]]))
    hp <- hp[is.finite(hp$private_bp), , drop = FALSE]
    hp$name   <- as.character(hp$haplotype)
    hp$is_ref <- hp$name == ref_hap
    nh  <- length(unique(hp$name))
    hgt <- max(3.8, 0.42 * nh + 2.4)
    fill_sc <- scale_fill_manual(values = c(`FALSE` = "#7fbf7b", `TRUE` = "#2c7fb8"),
                                 labels = c(`FALSE` = "haplotype", `TRUE` = "reference"),
                                 guide  = if (any(hp$is_ref)) "legend" else "none")
    # shared haplotype order from the full arm, so the two facets are comparable
    ordsrc <- if ("full" %in% hp$flavor) hp[hp$flavor == "full", ] else hp
    even   <- 100 / nh

    hp$y1 <- factor(hp$name, levels = ordsrc$name[order(ordsrc$pct_of_private)])
    p7 <- ggplot(hp, aes(pct_of_private, y1, fill = is_ref)) +
      geom_col(width = 0.72) + fill_sc +
      geom_vline(xintercept = even, linetype = 2, colour = "grey35") +
      geom_text(aes(label = mb(private_bp)), hjust = -0.12, size = 2.9, colour = "grey20") +
      facet_wrap(~ flavor, ncol = 2) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
      labs(title = paste0(label, " \u2014 who owns the private sequence"),
           subtitle = sprintf("dashed line = even share (%.1f%%). The reference is the graph backbone and is never clipped, so its RANK differs between arms.", even),
           x = "% of the cohort's private sequence", y = NULL, fill = NULL)
    ggsave(op(".private_by_haplotype.png"), p7, width = 11, height = hgt, dpi = 150)

    hp$y2 <- factor(hp$name, levels = ordsrc$name[order(ordsrc$pct_of_haplotype)])
    p8 <- ggplot(hp, aes(pct_of_haplotype, y2, fill = is_ref)) +
      geom_col(width = 0.72) + fill_sc +
      geom_text(aes(label = sprintf("%s / %s", mb(private_bp), mb(hap_bp))),
                hjust = -0.08, size = 2.7, colour = "grey20") +
      facet_wrap(~ flavor, ncol = 2) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.34))) +
      labs(title = paste0(label, " \u2014 private fraction of each haplotype"),
           subtitle = "private bp / this haplotype's total graph bp",
           x = "% of the haplotype that is private to it", y = NULL, fill = NULL)
    ggsave(op(".private_fraction_by_haplotype.png"), p8, width = 11, height = hgt, dpi = 150)
  }
}

# ---- audit ---------------------------------------------------------------------------
note("label", label)
note("min_bp", min_bp)
note("reference", ref_hap)
ad <- data.frame(metric = names(audit), value = unlist(audit, use.names = FALSE),
                 stringsAsFactors = FALSE)
write.table(ad[order(ad$metric), ], op(".private_figures_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[private_plots] %s: %d audit metrics, %d figures written",
                label, nrow(ad),
                length(list.files(outdir, pattern = paste0("^", label, "\\..*\\.png$")))))
