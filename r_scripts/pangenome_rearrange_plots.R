#!/usr/bin/env Rscript
# ======================================================================================
# pangenome_rearrange_plots.R
# Repo location: r_scripts/pangenome_rearrange_plots.R
#
# Layer 1 figures: inversions, duplications and rearrangement candidates from odgi untangle.
#
# WHY LAYER 1 IS THE PRIMARY REARRANGEMENT INSTRUMENT, NOT A SUPPLEMENT
# --------------------------------------------------------------------
# Inversions in this graph are overwhelmingly NOT bubble-representable: 29 path-explicit
# alleles and 411 alignment-rescued, against ~21.5 Mb of inverted sequence on chr10 alone from
# path projection. And `self.cov > 1` is the ONLY duplication signal available anywhere -- AT
# traversals found 27 node re-visits in 3,268,312 SV alleles. So these figures are not a second
# opinion on the variant catalog; for rearrangement they are the entire evidence base.
#
# WHAT THE CANDIDATE FIGURE HAS TO GET RIGHT
# ------------------------------------------
# Carriers are reported separately as chromosome-scale and unplaced, because an unplaced contig
# that projects inverted has no forward component to compare against and is therefore flagged
# BY CONSTRUCTION -- 74 of the 80 carriers at the chr10 locus were unplaced, and all 74 were
# flagged. Pooling them makes any_artifact_flag meaningless above a few carriers and buries the
# placed evidence, which is the part that distinguishes biology from assembly noise.
#
# n_chrom_individuals, not n_chrom_carriers, is the number that matters for interpretation:
# both haplotypes of one individual carrying a locus is ONE observation of that allele. The
# chr10 locus has six chromosome-scale carriers across only THREE individuals.
#
# Usage:
#   Rscript pangenome_rearrange_plots.R <label> <outdir> \
#       candidates=<tsv> inversions=<bed> duplications=<tsv> \
#       orientation=<tsv> audit=<tsv> [min_span=1000]
#
# INPUT SCHEMAS (verified against what rearrange_from_untangle.py writes)
#   candidates   ref_contig locus_start locus_end span_bp union_bp fill n_carriers
#                n_chrom_carriers n_unplaced_carriers n_chrom_individuals
#                chrom_carriers unplaced_carriers chrom_artifact_flag any_artifact_flag
#   inversions   BED: chrom start end name(query|inv) score strand   + a leading track line
#   duplications ref_contig ref_start ref_end query self_cov query_bp
#   orientation  query assembly contig fwd_bp inv_bp pct_inv pct_inv_unfiltered dup_bp
#                n_ref_targets harm_class harm_orient harm_flags artifact_flag
#   audit        metric value
#
# Outputs (in outdir):
#   <label>.rearrange_candidates.png     <label>.inverted_bp_by_chromosome.png
#   <label>.duplication_spectrum.png     <label>.orientation_flags.png
#   <label>.rearrange_figures_audit.tsv
# ======================================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: pangenome_rearrange_plots.R <label> <outdir> key=value ...")
label <- args[1]; outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
op <- function(x) file.path(outdir, paste0(label, x))
theme_set(theme_bw(base_size = 12))

kv <- list()
for (a in args[-(1:2)]) if (grepl("=", a, fixed = TRUE))
  kv[[sub("=.*$", "", a)]] <- sub("^[^=]*=", "", a)
optc <- function(k, d = "NONE") if (!is.null(kv[[k]]) && nzchar(kv[[k]])) kv[[k]] else d
optn <- function(k, d) { v <- suppressWarnings(as.numeric(optc(k, NA))); if (length(v) == 1 && is.finite(v)) v else d }

is_missing <- function(x) is.null(x) || is.na(x) || !nzchar(x) ||
  grepl("^(NONE|NO_)", basename(x)) || !file.exists(x) || file.size(x) == 0

# read.delim does NOT guarantee numeric types -- with check.names = FALSE and any stray
# non-numeric token a whole column comes back character, and cut() then fails with
# "'x' must be numeric". Coerce the columns we do arithmetic on, explicitly and on read,
# rather than trusting the type.
num_cols <- function(d, cols) {
  for (cc in intersect(cols, names(d))) d[[cc]] <- suppressWarnings(as.numeric(d[[cc]]))
  d
}

rd <- function(f, sep = "\t", header = TRUE, numeric_cols = character(0)) {
  if (is_missing(f)) return(NULL)
  d <- try(read.delim(f, sep = sep, header = header, comment.char = "#",
                      stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
  if (inherits(d, "try-error") || is.null(d) || !nrow(d)) return(NULL)
  if (length(numeric_cols)) d <- num_cols(d, numeric_cols)
  d
}

min_span <- optn("min_span", 1000)
audit <- list()
note  <- function(k, v) audit[[k]] <<- as.character(v)
skip  <- function(fig, why) { message(sprintf("[rearrange_plots] SKIP %s -- %s", fig, why))
                              note(paste0("skipped_", fig), why) }
mbp   <- function(x) ifelse(x >= 1e6, sprintf("%.1f Mb", x / 1e6), sprintf("%.0f kb", x / 1e3))

# chr10_1 and chr10_17+chr11_12 both pool to chr10, so pieces of one chromosome do not
# fragment an axis into dozens of near-empty categories
chrom_of <- function(x) sub("_.*$", "", sub("^.*#", "", x))
chrom_factor <- function(v) {
  u <- unique(v); n <- suppressWarnings(as.numeric(sub("^chr", "", u)))
  factor(v, levels = u[order(is.na(n), n, u)])
}

cand <- rd(optc("candidates"),
           numeric_cols = c("locus_start", "locus_end", "span_bp", "union_bp", "fill",
                            "n_carriers", "n_chrom_carriers", "n_unplaced_carriers",
                            "n_chrom_individuals"))
dups <- rd(optc("duplications"),
           numeric_cols = c("ref_start", "ref_end", "self_cov", "query_bp"))
ori  <- rd(optc("orientation"),
           numeric_cols = c("fwd_bp", "inv_bp", "pct_inv", "pct_inv_unfiltered",
                            "dup_bp", "n_ref_targets"))
aud  <- rd(optc("audit"))
inv  <- NULL
if (!is_missing(optc("inversions"))) {
  # BED has a leading `track name=...` line and no header row
  raw <- try(read.delim(optc("inversions"), sep = "\t", header = FALSE,
                        stringsAsFactors = FALSE, check.names = FALSE,
                        comment.char = ""), silent = TRUE)
  if (!inherits(raw, "try-error") && !is.null(raw) && nrow(raw)) {
    raw <- raw[!grepl("^track", raw[[1]]), , drop = FALSE]
    if (nrow(raw) && ncol(raw) >= 5) {
      inv <- data.frame(chrom = raw[[1]],
                        start = suppressWarnings(as.numeric(raw[[2]])),
                        end   = suppressWarnings(as.numeric(raw[[3]])),
                        name  = raw[[4]],
                        bp    = suppressWarnings(as.numeric(raw[[5]])),
                        stringsAsFactors = FALSE)
      inv <- inv[is.finite(inv$start) & is.finite(inv$end), , drop = FALSE]
      if (!nrow(inv)) inv <- NULL
    }
  }
}

# carry the audit's own numbers through, so the figure set and the table agree
if (!is.null(aud)) for (i in seq_len(nrow(aud)))
  note(paste0("untangle_", aud[[1]][i]), aud[[2]][i])

# ======================================================================================
# 1. Candidate loci: span vs chromosome-scale support
#
# The axis that matters is n_chrom_individuals, not n_carriers. 160 single-carrier loci exist
# and are almost all unplaced projections; the interesting locus is large AND carried by
# several individuals at chromosome scale. Point size shows how much unplaced support rides
# along, which is informative rather than disqualifying.
# ======================================================================================
if (is.null(cand)) {
  skip("rearrange_candidates", "no rearrangement_candidates.tsv")
} else if (!"n_chrom_individuals" %in% names(cand)) {
  skip("rearrange_candidates",
       "table predates the carrier split -- rerun REARRANGE with the current script")
} else {
  # a log10 y-axis silently discards span_bp <= 0, which produced a
  # "Removed 505 rows" warning on the first real run. Drop them explicitly and record the
  # count, so a large number is visible in the audit rather than buried in a warning.
  cd <- cand[is.finite(cand$span_bp) & cand$span_bp >= min_span, , drop = FALSE]
  note("candidate_rows_unplottable",
       sum(!is.finite(cand$span_bp) | cand$span_bp < min_span))
  if (!nrow(cd)) {
    skip("rearrange_candidates", sprintf("no loci >= %d bp", min_span))
  } else {
    cd$chrom  <- chrom_factor(chrom_of(cd$ref_contig))
    cd$flagged <- !is.na(cd$chrom_artifact_flag) &
                  cd$chrom_artifact_flag != "." & nzchar(cd$chrom_artifact_flag)
    cd$n_chrom_individuals[!is.finite(cd$n_chrom_individuals)] <- 0
    cd$n_unplaced_carriers[!is.finite(cd$n_unplaced_carriers)] <- 0
    note("candidate_loci", nrow(cd))
    note("candidate_loci_no_placed_carrier", sum(cd$n_chrom_carriers == 0))
    note("candidate_max_span_bp", max(cd$span_bp))

    p1 <- ggplot(cd, aes(n_chrom_individuals, span_bp)) +
      geom_point(aes(size = n_unplaced_carriers, colour = flagged), alpha = 0.75) +
      scale_y_log10(labels = function(v) mbp(v)) +
      scale_size_continuous(range = c(1.2, 6)) +
      scale_colour_manual(values = c(`FALSE` = "#2c7fb8", `TRUE` = "firebrick"),
                          labels = c(`FALSE` = "clean", `TRUE` = "a placed carrier flagged")) +
      labs(title = paste0(label, " — rearrangement candidate loci"),
           subtitle = paste0("x = INDIVIDUALS with chromosome-scale support (both haplotypes ",
                             "of one individual = one observation)"),
           x = "individuals with chromosome-scale carriers", y = "locus span (log)",
           size = "unplaced\ncarriers", colour = NULL) +
      theme(legend.position = "right")
    ggsave(op(".rearrange_candidates.png"), p1, width = 9, height = 5.5, dpi = 150)
  }
}

# ======================================================================================
# 2. Inverted bp per reference chromosome
#
# Where the inverted sequence actually is. chr10 dominates on this cohort, which is the
# observation the candidate table exists to explain.
# ======================================================================================
if (is.null(inv)) {
  skip("inverted_bp_by_chromosome", "no inversions.bed rows")
} else {
  iv <- inv
  iv$len <- iv$end - iv$start
  iv <- iv[is.finite(iv$len) & iv$len >= min_span, , drop = FALSE]
  if (!nrow(iv)) {
    skip("inverted_bp_by_chromosome", sprintf("no segments >= %d bp", min_span))
  } else {
    # name is "<query>|inv"; the query's own assembly is the carrier
    iv$carrier <- sub("\\|inv$", "", iv$name)
    iv$carrier <- sub("#.*$", "", iv$carrier)
    iv$chrom <- chrom_factor(chrom_of(iv$chrom))
    agg <- aggregate(len ~ chrom + carrier, data = iv, FUN = sum)
    note("inverted_bp_total", sum(iv$len))
    note("inverted_segments", nrow(iv))

    p2 <- ggplot(agg, aes(chrom, len, fill = carrier)) +
      geom_col() +
      scale_y_continuous(labels = function(v) mbp(v)) +
      labs(title = paste0(label, " — inverted sequence per reference chromosome"),
           subtitle = paste0("odgi untangle minus-strand projections >= ", min_span,
                             " bp, stacked by carrying assembly"),
           x = NULL, y = "inverted sequence", fill = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right", legend.text = element_text(size = 8))
    ggsave(op(".inverted_bp_by_chromosome.png"), p2, width = 10, height = 5.5, dpi = 150)
  }
}

# ======================================================================================
# 3. Duplication size spectrum
#
# self.cov > 1 is the ONLY duplication signal in the pipeline. AT traversals found 27 node
# re-visits in 3.27M SV alleles, so nothing else can produce this figure.
# ======================================================================================
if (is.null(dups)) {
  skip("duplication_spectrum", "no duplications.tsv")
} else {
  dp <- dups[is.finite(dups$query_bp) & dups$query_bp >= min_span &
             is.finite(dups$self_cov), , drop = FALSE]
  note("duplication_rows_unplottable",
       sum(!is.finite(dups$query_bp) | !is.finite(dups$self_cov)))
  if (!nrow(dp)) {
    skip("duplication_spectrum", sprintf("no duplications >= %d bp", min_span))
  } else {
    dp$chrom <- chrom_factor(chrom_of(dp$ref_contig))
    dp$cov_bin <- cut(dp$self_cov, breaks = c(1, 2, 3, 5, 10, Inf),
                      labels = c("2", "3", "4-5", "6-10", ">10"),
                      right = TRUE, include.lowest = FALSE)
    note("duplication_loci", nrow(dp))
    note("duplication_bp_total", sum(dp$query_bp))
    note("duplication_max_self_cov", max(dp$self_cov))

    p3 <- ggplot(dp, aes(query_bp, fill = cov_bin)) +
      geom_histogram(bins = 40) +
      scale_x_log10(labels = function(v) mbp(v)) +
      labs(title = paste0(label, " — duplication size spectrum"),
           subtitle = paste0("self.cov > 1 from untangle: the ONLY duplication signal ",
                             "available (AT found 27 re-visits in 3.27M alleles)"),
           x = "duplicated query span (log)", y = "loci", fill = "self.cov") +
      theme(legend.position = "right")
    ggsave(op(".duplication_spectrum.png"), p3, width = 9, height = 5, dpi = 150)
  }
}

# ======================================================================================
# 4. Orientation flags by query class
#
# ORIENTATION_SUSPECT fires on 2,615 queries, nearly all unplaced. This figure exists to make
# that visible so the flag is not read as 2,615 orientation errors: an unplaced contig that
# projects inverted has no forward component, so pct_inv is ~100% by construction.
#
# pct_inv vs pct_inv_unfiltered is the diagnostic -- a large gap means the -j / size filters
# are driving the number rather than orientation, which is why ORIENTATION_SUSPECT requires
# BOTH to be high.
# ======================================================================================
# `o` must exist before the branch: the guard below references it, and if `ori` is NULL an
# unassigned `o` is an "object not found" error rather than a skipped figure.
o <- NULL
if (is.null(ori)) {
  skip("orientation_flags", "no query_orientation.tsv")
} else {
  o <- ori
  o <- o[is.finite(o$pct_inv) & is.finite(o$pct_inv_unfiltered), , drop = FALSE]
  if (!nrow(o)) {
    skip("orientation_flags", "no rows with finite pct_inv")
    o <- NULL
  }
}
if (!is.null(o) && nrow(o)) {
  o$placed <- ifelse(grepl("^unplaced", o$contig), "unplaced", "placed")
  o$flag <- ifelse(is.na(o$artifact_flag) | o$artifact_flag == "." | !nzchar(o$artifact_flag),
                   "none", o$artifact_flag)
  note("orientation_queries", nrow(o))
  note("orientation_flagged", sum(o$flag != "none"))
  note("orientation_flagged_unplaced",
       sum(o$flag != "none" & o$placed == "unplaced"))

  p4a <- ggplot(o, aes(pct_inv_unfiltered, pct_inv, colour = placed)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_point(alpha = 0.5, size = 1) +
    scale_colour_manual(values = c(placed = "#2c7fb8", unplaced = "grey50")) +
    labs(title = paste0(label, " — orientation: filtered vs unfiltered"),
         subtitle = paste0("off-diagonal means the -j / size filters drive the number, ",
                           "not orientation. ORIENTATION_SUSPECT needs BOTH high."),
         x = "% inverted, all rows", y = "% inverted, filtered rows", colour = NULL) +
    theme(legend.position = "top")
  ggsave(op(".orientation_flags.png"), p4a, width = 8, height = 5.5, dpi = 150)

  ft <- as.data.frame(table(flag = o$flag, placed = o$placed), stringsAsFactors = FALSE)
  ft <- ft[ft$Freq > 0, , drop = FALSE]
  if (nrow(ft)) {
    p4b <- ggplot(ft, aes(reorder(flag, Freq), Freq, fill = placed)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c(placed = "#2c7fb8", unplaced = "grey50")) +
      labs(title = paste0(label, " — artifact flags by query class"),
           subtitle = paste0("an unplaced contig projecting inverted has no forward ",
                             "component, so it is flagged by construction"),
           x = NULL, y = "queries", fill = NULL) +
      theme(legend.position = "top")
    ggsave(op(".orientation_flag_counts.png"), p4b, width = 8,
           height = max(3.5, 0.4 * nrow(ft) + 2.4), dpi = 150)
  }
}

# ---- audit ---------------------------------------------------------------------------
note("label", label)
note("min_span", min_span)
ad <- data.frame(metric = names(audit), value = unlist(audit, use.names = FALSE),
                 stringsAsFactors = FALSE)
write.table(ad[order(ad$metric), ], op(".rearrange_figures_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[rearrange_plots] %s: %d audit metrics, %d figures written",
                label, nrow(ad),
                length(list.files(outdir, pattern = paste0("^", label, "\\..*\\.png$")))))
