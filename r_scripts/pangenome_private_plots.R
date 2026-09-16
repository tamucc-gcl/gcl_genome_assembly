#!/usr/bin/env Rscript
# ======================================================================================
# pangenome_private_plots.R
# Repo location: r_scripts/pangenome_private_plots.R
#
# The private-sequence figure set, answering Chris Bird's requests (2026-08-26, 2026-09):
#
#   "for the SV histogram, it would be interesting to make an alternative version with
#    coarse bins on the x axis and the y axis representing number of bp in SV, rather than
#    number of SV"
#   "a similar plot for private haplotypes would also be informative"
#   "It might be a good idea to independently map the private haplotypes against the other
#    assemblies"
#   "The private haplotypes could also be driven by transposons"
#
# EVERYTHING IS bp-WEIGHTED. That is the whole point of the request, and it matters: 98.9%
# of private segments are under 100 bp and carry 4% of the private sequence, so a
# count-weighted figure is a spike at the left edge that answers a different question.
#
# Usage:
#   Rscript pangenome_private_plots.R <label> <outdir> \
#       spectrum_clip= spectrum_full= tier_clip= sv_clip= \
#       hap_clip= hap_full= contig_clip= contig_full= \
#       evidence_clip= evidence_full= xtab_clip= xtab_full= \
#       [ref=] [bins=] [min_bp=1000]
#
# Every input is optional: a missing table skips its figure with a recorded reason rather
# than failing the task.
#
# INPUT SCHEMAS, verified against real output
#   private_segment_spectrum.tsv    scope size_bin n_segments segment_bp
#   coverage_segment_spectrum.tsv   scope tier size_bin n_segments segment_bp
#   size_spectrum.tsv               primary_class size_bin n_alleles per_allele_bp
#   hap_private.tsv                 haplotype sample hap private_bp hap_bp
#                                   pct_of_private pct_of_haplotype
#   hap_private_by_contig.tsv       haplotype contig private_bp contig_graph_bp
#                                   pct_of_contig pct_of_hap_private pct_of_pangenome_private
#   private_evidence.csv            27 columns; see PANGENOME_PRIVATE_JOIN
#   private_evidence_xtab.tsv       scope key set combined n_segments segment_bp
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
optn <- function(k, d) { v <- suppressWarnings(as.numeric(optc(k, NA)))
                         if (length(v) == 1 && is.finite(v)) v else d }

is_missing <- function(x) is.null(x) || is.na(x) || !nzchar(x) ||
  grepl("^(NONE|NO_)", basename(x)) || !file.exists(x) || file.size(x) == 0

# Comment lines are counted and skipped, NOT handled with comment.char. '#' is a legitimate
# character inside a PanSN haplotype name (Sde-CBau_104#1), and comment.char = "#" truncates
# every key at the separator and shifts every subsequent column left -- which made a
# per-chromosome figure report "no chr* contigs" from a file holding 602 of them.
read_nohash <- function(f, sep = "\t") {
  n_skip <- 0L
  con <- file(f, "r")
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (!length(ln)) break
    if (!grepl("^[[:space:]]*#", ln)) break
    n_skip <- n_skip + 1L
  }
  close(con)
  d <- try(read.delim(f, sep = sep, header = TRUE, skip = n_skip,
                      stringsAsFactors = FALSE, check.names = FALSE,
                      comment.char = "", quote = ""), silent = TRUE)
  if (inherits(d, "try-error") || is.null(d) || !nrow(d)) return(NULL)
  d
}

# read.delim does not guarantee numeric types: one stray token turns a column character and
# the failure surfaces at plot time, not read time.
num_cols <- function(d, cols) {
  for (cc in intersect(cols, names(d))) d[[cc]] <- suppressWarnings(as.numeric(d[[cc]]))
  d
}
rd <- function(f, sep = "\t", numeric_cols = character(0)) {
  if (is_missing(f)) return(NULL)
  d <- read_nohash(f, sep = sep)
  if (is.null(d)) return(NULL)
  if (length(numeric_cols)) d <- num_cols(d, numeric_cols)
  d
}

ref_hap <- optc("ref", "NONE")
min_bp  <- optn("min_bp", 1000)
same_hap <- function(key, ref) {
  if (is.null(ref) || ref == "NONE") return(rep(FALSE, length(key)))
  sub("#.*$", "", as.character(key)) == sub("#.*$", "", ref) | as.character(key) == ref
}

gbp <- function(x) ifelse(is.na(x), "NA",
                   ifelse(x >= 1e9, sprintf("%.2f Gb", x / 1e9),
                   ifelse(x >= 1e6, sprintf("%.0f Mb", x / 1e6),
                          sprintf("%.0f kb", x / 1e3))))
mb <- function(x) sprintf("%.1f Mb", x / 1e6)

audit <- list()
note  <- function(k, v) audit[[k]] <<- as.character(v)
skip  <- function(fig, why) { message(sprintf("[private_plots] SKIP %s -- %s", fig, why))
                              note(paste0("skipped_", fig), why) }
has_facet <- function(d, col) !is.null(d) && nrow(d) > 0 && col %in% names(d) &&
  length(unique(stats::na.omit(d[[col]]))) > 0

# size_bin is a string like "1000-5000" or ">=1000000"; order by its lower edge or 100-500
# sorts before 0-100
bin_order <- function(x) {
  lo <- suppressWarnings(as.numeric(sub("^>=", "", sub("-.*$", "", x))))
  factor(x, levels = unique(x[order(lo)]))
}

bind_flavours <- function(clip_f, full_f, reader) {
  parts <- list()
  for (fl in c("clip", "full")) {
    d <- reader(if (fl == "clip") clip_f else full_f)
    if (!is.null(d)) { d$flavor <- fl; parts[[fl]] <- d }
  }
  if (!length(parts)) return(NULL)
  do.call(rbind, parts)
}

spec <- bind_flavours(optc("spectrum_clip"), optc("spectrum_full"),
        function(f) rd(f, numeric_cols = c("n_segments", "segment_bp")))
hap  <- bind_flavours(optc("hap_clip"), optc("hap_full"),
        function(f) rd(f, numeric_cols = c("private_bp", "hap_bp",
                                "pct_of_private", "pct_of_haplotype")))
ctg  <- bind_flavours(optc("contig_clip"), optc("contig_full"),
        function(f) rd(f, numeric_cols = c("private_bp", "contig_graph_bp",
                                "pct_of_contig")))
xtab <- bind_flavours(optc("xtab_clip"), optc("xtab_full"),
        function(f) rd(f, numeric_cols = c("n_segments", "segment_bp")))
ev   <- bind_flavours(optc("evidence_clip"), optc("evidence_full"),
        function(f) rd(f, sep = ",",
                       numeric_cols = c("span_bp", "is_private", "repeat_like",
                                        "median_copy", "copy_ratio",
                                        "aligned_frac_merged", "n_other_assemblies")))
tier <- rd(optc("tier_clip"), numeric_cols = c("n_segments", "segment_bp"))
tcon <- rd(optc("tier_contig_clip"),
           numeric_cols = c("bp", "contig_total_bp", "pct_of_contig"))
svsp <- rd(optc("sv_clip"),   numeric_cols = c("n_alleles", "per_allele_bp"))

# ======================================================================================
# FIGURE 1 (was 2). Where the sequence is, by segment size AND by how shared it is.
#
# THIS IS THE DIRECT ANSWER to "can we make a private haplotype size histogram? I'd expect
# it to be similar to the SV". Two panels on one x axis:
#
#   top     SV bp per size bin, black -- the histogram that already existed, re-weighted
#           from count to bp as asked
#   bottom  graph bp per size bin, split by SHARING TIER (core / soft-core / shell /
#           private), so the private distribution can be read against the shared ones
#
# CLIP ONLY, deliberately. The SV catalog exists only on the clip graph, so a full-arm
# bottom panel would not be comparable to the top one -- and the cactus authors recommend
# the clip graph for most applications. The full arm's extra private sequence is also
# inflated by an artifact: a 111 Mb scaffold in Sde-CTlk_104_hap1 fuses chr5 and chr9, and
# cactus assigns each contig to ONE chromosome, so ~73 Mb of unaligned chr9 sits in chr5's
# subgraph at coverage 1. Clipping removes exactly that class of sequence.
# ======================================================================================
if (is.null(tier)) {
  skip("size_by_sharing", "no coverage_segment_spectrum.tsv for the clip graph")
} else {
  tt <- tier[tier$scope == "ALL" & is.finite(tier$segment_bp), , drop = FALSE]
  if (!nrow(tt)) {
    skip("size_by_sharing", "no ALL-scope rows in the tier spectrum")
  } else {
    lv <- c("core", "soft-core", "shell", "private")
    tt$tier <- factor(tt$tier, levels = lv[lv %in% unique(tt$tier)])
    allb <- unique(c(as.character(tt$size_bin),
                     if (!is.null(svsp)) as.character(svsp$size_bin) else character(0)))
    allb <- allb[allb != "0-50"]
    ord <- levels(bin_order(allb))
    tt$size_bin <- factor(as.character(tt$size_bin), levels = ord)

    for (l in levels(tt$tier))
      note(paste0("tier_bp.", l), sum(tt$segment_bp[tt$tier == l], na.rm = TRUE))

    # DODGED, not stacked. Stacked, a reader cannot compare tiers within a bin -- the whole
    # question is whether private sequence sits at different sizes than shared sequence, and
    # that is a comparison between bars, not a total.
    # BOTH PANELS PROPORTIONAL, and both floored at 50 bp. They were not comparable before:
    # the SV panel is per-allele bp totalling 2.97 Gb with a 50 bp floor (SNPs and indels
    # excluded), the tier panel deduplicated graph bp totalling 1.93 Gb with no floor. Chris
    # asked whether the SHAPES match, so share of own total is the right axis and the ranges
    # have to agree.
    tt <- tt[!(as.character(tt$size_bin) %in% c("0-50")), , drop = FALSE]
    tt$frac <- tt$segment_bp / sum(tt$segment_bp)

    p_bot <- ggplot(tt, aes(size_bin, frac, fill = tier)) +
      geom_col(position = position_dodge(preserve = "single")) +
      scale_y_continuous(labels = function(v) sprintf("%.0f%%", 100 * v)) +
      scale_x_discrete(limits = ord, drop = FALSE) +
      scale_fill_manual(values = c(core = "#2166ac", `soft-core` = "#67a9cf",
                                   shell = "#fdae61", private = "#d73027")) +
      labs(x = "segment size (bp)", y = "share of graph bp", fill = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

    if (is.null(svsp)) {
      skip("size_by_sharing_sv_panel", "no SV size_spectrum.tsv; bottom panel only")
      ggsave(op(".size_by_sharing.png"), p_bot +
               labs(title = paste0(label, " \u2014 where the sequence is")),
             width = 9, height = 5.5, dpi = 150)
    } else {
      # SV classes only. SNP and INDEL would dominate the count but are not what
      # "sequence in structural variants" means.
      # Guard the schema. This wants the PRE-AGGREGATED size_spectrum.tsv
      # (primary_class / size_bin / n_alleles / per_allele_bp), NOT sv_sizes.tsv, which is
      # one row per allele -- 31.3M of them. Passing the wrong one gave
      # "object 'per_allele_bp' not found" from inside aggregate(), which says nothing about
      # which file was wrong.
      need <- c("primary_class", "size_bin", "per_allele_bp")
      if (!all(need %in% names(svsp))) {
        skip("size_by_sharing_sv_panel",
             sprintf("sv_clip has columns [%s]; expected size_spectrum.tsv with [%s]",
                     paste(names(svsp), collapse = ","), paste(need, collapse = ",")))
        svsp <- NULL
      }
      sv <- if (is.null(svsp)) NULL else
        svsp[svsp$primary_class %in%
               c("INS", "DEL", "SUBST", "DUP", "INV_DUP", "INV_PATH_EXPLICIT") &
             is.finite(svsp$per_allele_bp), , drop = FALSE]
      if (is.null(sv) || !nrow(sv)) {
        ggsave(op(".size_by_sharing.png"), p_bot +
                 labs(title = paste0(label, " \u2014 where the sequence is")),
               width = 9, height = 5.5, dpi = 150)
        sv <- NULL
      }
      agg <- if (is.null(sv)) NULL else aggregate(per_allele_bp ~ size_bin, data = sv, FUN = sum)
      if (!is.null(agg)) {
      agg <- agg[!(as.character(agg$size_bin) %in% c("0-50")), , drop = FALSE]
      agg$frac <- agg$per_allele_bp / sum(agg$per_allele_bp)
      agg$size_bin <- factor(as.character(agg$size_bin), levels = ord)
      note("sv_bp_total", sum(agg$per_allele_bp))
      # per_allele_bp sums EVERY alt allele's length, and SUBST is ~60% of the catalog with
      # up to seven alts at a locus -- so this totals 2.97 Gb against a ~1 Gb reference. The
      # SHAPE across bins is meaningful; the total is not a genome quantity, unlike the
      # deduplicated graph bp below it. Recorded so the asymmetry is not read as a result.
      note("sv_bp_note", "per-allele sum; double-counts multi-allelic loci, unlike the graph bp below")

      p_top <- ggplot(agg, aes(size_bin, frac)) +
        geom_col(fill = "grey15") +
        scale_y_continuous(labels = function(v) sprintf("%.0f%%", 100 * v)) +
        scale_x_discrete(limits = ord, drop = FALSE) +
        labs(title = paste0(label, " \u2014 size distribution: SVs vs graph sharing"),
             subtitle = "share of each panel's own total, >=50 bp, clip graph",
             x = NULL, y = "share of SV bp") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

      # stacked without a layout dependency: write both panels and a combined image via
      # gridExtra only if it is available, else two files
      okg <- requireNamespace("gridExtra", quietly = TRUE)
      if (okg) {
        g <- gridExtra::arrangeGrob(p_top, p_bot, ncol = 1, heights = c(1, 1.6))
        ggsave(op(".size_by_sharing.png"), g, width = 9, height = 8.5, dpi = 150)
      } else {
        note("gridExtra", "absent -- panels written separately")
        ggsave(op(".size_by_sharing_sv.png"), p_top, width = 9, height = 4, dpi = 150)
        ggsave(op(".size_by_sharing.png"),
               p_bot + labs(title = paste0(label, " \u2014 where the sequence is")),
               width = 9, height = 5.5, dpi = 150)
      }
      }
    }
  }
}

# ======================================================================================
# FIGURE 2. Private vs its matched control, in bp.
#
# bp-WEIGHTED, not segment-weighted. The two agree closely here (65.0% of segments vs
# 62.6% of bp for PRIVATE_CONFIRMED+REPEAT_LIKE), and that agreement is itself worth
# showing -- it says the result is not an artifact of counting many tiny segments.
#
# The control is size-matched, cross-individual, NON-private windows from the same
# haplotype, run through the identical pipeline. Without it neither measure means anything:
# "private" as the graph defines it (a node walked by one haplotype) is not the claim that
# the sequence exists nowhere else.
# ======================================================================================
if (!has_facet(ev, "flavor") || !has_facet(ev, "set")) {
  skip("private_vs_control", "no private_evidence.csv rows with flavour and set")
} else {
  rows <- list()
  for (fl in unique(ev$flavor)) for (st in unique(ev$set)) {
    s <- ev[ev$flavor == fl & ev$set == st, , drop = FALSE]
    if (!nrow(s)) next
    tb <- sum(s$span_bp, na.rm = TRUE)
    if (!is.finite(tb) || tb <= 0) next
    rl <- is.finite(s$repeat_like) & s$repeat_like == 1
    np <- s$map_verdict == "NOT_PRIVATE"
    rows[[length(rows) + 1]] <- data.frame(
      flavor = fl, set = st, bp = tb, measure = "repeat-like (k-mer copy number)",
      value = sum(s$span_bp[rl], na.rm = TRUE) / tb, stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      flavor = fl, set = st, bp = tb, measure = "aligns to another assembly",
      value = sum(s$span_bp[np], na.rm = TRUE) / tb, stringsAsFactors = FALSE)
  }
  cmp <- do.call(rbind, rows)
  if (is.null(cmp) || !nrow(cmp)) {
    skip("private_vs_control", "no private/control rows with usable bp")
  } else {
    p2 <- ggplot(cmp, aes(set, value, fill = set)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = sprintf("%.1f%%\n(%s)", 100 * value, gbp(bp))),
                vjust = -0.15, size = 3) +
      facet_grid(measure ~ flavor) +
      scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.25),
                         labels = function(v) sprintf("%.0f%%", 100 * v)) +
      scale_fill_manual(values = c(private = "#d73027", control = "grey60")) +
      labs(title = paste0(label, " \u2014 private sequence vs a matched control"),
           subtitle = "share of bp; control = matched non-private windows",
           x = NULL, y = "share of bp") +
      theme(legend.position = "none")
    ggsave(op(".private_vs_control.png"), p2, width = 8, height = 6, dpi = 150)
    for (i in seq_len(nrow(cmp)))
      note(sprintf("%s.%s.%s", gsub("[^a-z]", "", tolower(cmp$measure[i])),
                   cmp$flavor[i], cmp$set[i]), sprintf("%.4f", cmp$value[i]))
  }
}

# ======================================================================================
# FIGURE 3. Copy-number distribution, private vs control.
#
# THE STRONGEST EVIDENCE FOR THE TRANSPOSON ANSWER, and better than any single threshold
# statistic. Measured on the clip graph:
#
#   copy ratio    control   private
#   0.5-1.5        43.7%      2.8%     <- control peaks at SINGLE COPY
#   10-100         11.8%     45.4%
#   >=100           3.9%     31.6%
#
# The control peaking at 1.0 is the calibration verified independently: the single-copy
# reference is derived from the control itself, so if it were set too low the control would
# shift with it. It does not. And private sequence sits at 10-100x and beyond -- satellite
# and high-copy TE territory, which is exactly the sequence that fails to align between
# haplotypes. Two distributions that barely overlap.
#
# copy_ratio = the segment's median k-mer multiplicity in its own sample's reads, divided
# by that haplotype's single-copy level. 1.0 means present once in the genome.
# ======================================================================================
if (is.null(ev) || !"copy_ratio" %in% names(ev)) {
  skip("copy_ratio", "no copy_ratio column in the evidence table")
} else {
  cr <- ev[is.finite(ev$copy_ratio) & ev$copy_ratio > 0, , drop = FALSE]
  if (!nrow(cr)) {
    skip("copy_ratio", "no rows with a usable copy_ratio")
  } else {
    brk <- c(0, 0.5, 1.5, 3, 10, 100, Inf)
    lab <- c("<0.5", "0.5-1.5", "1.5-3", "3-10", "10-100", ">=100")
    cr$bin <- cut(cr$copy_ratio, breaks = brk, labels = lab, right = FALSE)
    agg <- aggregate(span_bp ~ bin + set + flavor, data = cr, FUN = sum)
    tot <- aggregate(span_bp ~ set + flavor, data = cr, FUN = sum)
    names(tot)[names(tot) == "span_bp"] <- "tot"
    agg <- merge(agg, tot, by = c("set", "flavor"))
    agg$frac <- agg$span_bp / agg$tot

    p3 <- ggplot(agg, aes(bin, frac, fill = set)) +
      geom_col(position = position_dodge(preserve = "single")) +
      facet_wrap(~ flavor, ncol = 1) +
      scale_y_continuous(labels = function(v) sprintf("%.0f%%", 100 * v)) +
      scale_fill_manual(values = c(private = "#d73027", control = "grey60")) +
      labs(title = paste0(label, " \u2014 k-mer copy number, private vs control"),
           subtitle = "copies per haploid genome; 2.0 = normal diploid",
           x = "copy ratio (1.0 = single copy)", y = "share of bp", fill = NULL) +
      theme(legend.position = "top")
    ggsave(op(".copy_ratio.png"), p3, width = 9, height = 6, dpi = 150)

    for (i in seq_len(nrow(agg)))
      if (agg$flavor[i] == "clip")
        note(sprintf("copyratio.%s.%s", agg$set[i], gsub("[^0-9a-z.]", "", agg$bin[i])),
             sprintf("%.4f", agg$frac[i]))
  }
}

# ======================================================================================
# FIGURE 4. Evidence cross-tabulation per chromosome, in ABSOLUTE bp.
#
# NOT normalised to 100% per chromosome. The normalised version could not be read for how
# much private sequence a chromosome HAS -- a chromosome 5% private and one 25% private
# looked identical, and the 65% blue band invited reading as "65% of the chromosome" when
# it is 65% of the private subset. Absolute bp makes both the composition and the amount
# visible on one axis.
#
# The four combinations each mean something:
#   NOT_PRIVATE + REPEAT_LIKE        present elsewhere AND high copy -> graph collapse
#   PRIVATE_CONFIRMED + UNIQUE_LIKE  absent elsewhere AND single copy -> novel sequence
#   NOT_PRIVATE + UNIQUE_LIKE        the graph failed to merge homologous sequence
#   PRIVATE_CONFIRMED + REPEAT_LIKE  haplotype-specific high-copy array
# ======================================================================================
if (!has_facet(xtab, "flavor")) {
  skip("evidence_by_chromosome", "no private_evidence_xtab.tsv rows")
} else {
  x <- xtab[xtab$scope == "CHROM" & xtab$set == "private" &
            is.finite(xtab$segment_bp) & xtab$segment_bp > 0, , drop = FALSE]
  if (!nrow(x)) {
    skip("evidence_by_chromosome",
         "no CHROM-scope private rows with bp -- the xtab writes bp only for ALL scope")
  } else {
    ord <- unique(x$key)
    num <- suppressWarnings(as.numeric(sub("^chr", "", ord)))
    x$key <- factor(x$key, levels = ord[order(is.na(num), num, ord)])
    p4 <- ggplot(x, aes(key, segment_bp, fill = combined)) +
      geom_col() +
      facet_wrap(~ flavor, ncol = 1, scales = "free_y") +
      scale_y_continuous(labels = gbp) +
      labs(title = paste0(label, " \u2014 private sequence per chromosome, by evidence"),
           subtitle = "bp summed over all 10 haplotypes",
           x = NULL, y = "private sequence, summed over haplotypes", fill = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right", legend.text = element_text(size = 8))
    ggsave(op(".evidence_by_chromosome.png"), p4, width = 10, height = 7, dpi = 150)
  }
}

# ======================================================================================
# FIGURE 4b. Tier composition per chromosome -- private as a FRACTION of the chromosome.
#
# The bar height IS the chromosome's graph bp, split by sharing level, so a private share is
# readable as a slice of the chromosome. The per-chromosome EVIDENCE figure above sums private
# bp across all ten haplotypes, which made chr1 look like ~60 Mb of private sequence on a
# 94.7 Mb chromosome -- i.e. as though most of the chromosome sat in one haplotype. It is
# ~6 Mb per haplotype, about 6%. Same data, and only this framing cannot be misread.
# ======================================================================================
if (is.null(tcon)) {
  skip("tier_by_chromosome", "no coverage_by_contig.tsv for the clip graph")
} else {
  tc <- tcon[grepl("^chr", tcon$contig) & is.finite(tcon$bp), , drop = FALSE]
  if (!nrow(tc)) {
    skip("tier_by_chromosome", "no chr* contigs")
  } else {
    # chr10_1 and chr10_17+chr11_12 both pool to chr10
    tc$chrom <- sub("_.*$", "", tc$contig)
    agg <- aggregate(bp ~ chrom + tier, data = tc, FUN = sum)
    ordc <- unique(agg$chrom); nn <- suppressWarnings(as.numeric(sub("^chr", "", ordc)))
    agg$chrom <- factor(agg$chrom, levels = ordc[order(is.na(nn), nn, ordc)])
    lv <- c("core", "soft-core", "shell", "private")
    agg$tier <- factor(agg$tier, levels = lv[lv %in% unique(agg$tier)])

    tot <- aggregate(bp ~ chrom, data = agg, FUN = sum)
    pv <- agg[agg$tier == "private", c("chrom", "bp")]
    names(pv)[2] <- "priv"
    pv <- merge(tot, pv, by = "chrom", all.x = TRUE)
    pv$priv[is.na(pv$priv)] <- 0
    note("private_pct_of_chrom_range",
         sprintf("%.1f%%-%.1f%%", 100 * min(pv$priv / pv$bp), 100 * max(pv$priv / pv$bp)))

    p4b <- ggplot(agg, aes(chrom, bp, fill = tier)) +
      geom_col() +
      scale_y_continuous(labels = gbp) +
      scale_fill_manual(values = c(core = "#2166ac", `soft-core` = "#67a9cf",
                                   shell = "#fdae61", private = "#d73027")) +
      labs(title = paste0(label, " \u2014 chromosome composition by sharing level"),
           subtitle = "bar height is the chromosome; each node counted once",
           x = NULL, y = "graph sequence", fill = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
    ggsave(op(".tier_by_chromosome.png"), p4b, width = 10, height = 5.5, dpi = 150)
  }
}

# ======================================================================================
# FIGURE 5. Private fraction per chromosome.
#
# pct_of_contig, so it is a FRACTION question rather than a chromosome-length question.
# ======================================================================================
if (is.null(ctg)) {
  skip("private_by_chromosome", "no hap_private_by_contig.tsv")
} else {
  c2 <- ctg[grepl("^chr", ctg$contig) & is.finite(ctg$pct_of_contig), , drop = FALSE]
  if (!nrow(c2)) {
    skip("private_by_chromosome", "no chr* contigs")
  } else {
    c2$chrom <- sub("_.*$", "", c2$contig)
    ord <- unique(c2$chrom); num <- suppressWarnings(as.numeric(sub("^chr", "", ord)))
    c2$chrom <- factor(c2$chrom, levels = ord[order(is.na(num), num, ord)])
    p5 <- ggplot(c2, aes(chrom, pct_of_contig, colour = flavor)) +
      geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.75)) +
      scale_colour_manual(values = c(clip = "grey40", full = "#2c7fb8")) +
      labs(title = paste0(label, " \u2014 private fraction per chromosome"),
           subtitle = "one point per haplotype",
           x = NULL, y = "private share of the chromosome (%)", colour = "graph") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
    ggsave(op(".private_by_chromosome.png"), p5, width = 9, height = 5, dpi = 150)
  }
}

# ======================================================================================
# FIGURE 6-7. Ownership, and the clip-vs-full disagreement.
# ======================================================================================
if (!has_facet(hap, "flavor")) {
  skip("private_ownership", "no hap_private.tsv rows")
} else {
  h <- hap[is.finite(hap$private_bp), , drop = FALSE]
  h$is_ref <- same_hap(h$haplotype, ref_hap)
  o <- if ("full" %in% h$flavor) h[h$flavor == "full", ] else h
  h$haplotype <- factor(h$haplotype, levels = o$haplotype[order(o$private_bp)])
  nh <- length(unique(as.character(h$haplotype)))

  p6 <- ggplot(h, aes(private_bp, haplotype, fill = is_ref)) +
    geom_col(position = position_dodge(preserve = "single")) +
    facet_wrap(~ flavor, ncol = 2) +
    scale_x_continuous(labels = gbp) +
    scale_fill_manual(values = c(`FALSE` = "#7fbf7b", `TRUE` = "#2c7fb8"),
                      labels = c(`FALSE` = "haplotype", `TRUE` = "reference"),
                      guide = if (any(h$is_ref)) "legend" else "none") +
    labs(title = paste0(label, " \u2014 private sequence per haplotype"),
         subtitle = "the reference is never clipped, so its rank differs by graph",
         x = "private sequence", y = NULL, fill = NULL) +
    theme(legend.position = "top")
  ggsave(op(".private_clip_vs_full.png"), p6, width = 11,
         height = max(4, 0.42 * nh + 2.4), dpi = 150)

  for (fl in unique(h$flavor)) {
    hh <- h[h$flavor == fl, ]
    if (nrow(hh) && ref_hap != "NONE" && any(same_hap(hh$haplotype, ref_hap)))
      note(paste0("reference_rank_", fl),
           sprintf("%d of %d", rank(-hh$private_bp)[same_hap(hh$haplotype, ref_hap)], nrow(hh)))
  }
}

# ---- audit ---------------------------------------------------------------------------
note("label", label); note("min_bp", min_bp); note("reference", ref_hap)
ad <- data.frame(metric = names(audit), value = unlist(audit, use.names = FALSE),
                 stringsAsFactors = FALSE)
write.table(ad[order(ad$metric), ], op(".private_figures_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[private_plots] %s: %d metrics, %d figures",
                label, nrow(ad),
                length(list.files(outdir, pattern = paste0("^", label, "\\..*\\.png$")))))
