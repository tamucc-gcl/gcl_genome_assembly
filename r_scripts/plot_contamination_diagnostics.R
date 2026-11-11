library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
evidence <- read.delim(args[1], check.names = FALSE)
dir.create(args[2], showWarnings = FALSE)
evidence$gc_pct <- 100 * evidence$gc_fraction
special <- c("no-hit" = "#B1B7BD", "unresolved" = "#666E76", "conflicting" = "#A03B83")
taxon_span <- sort(tapply(evidence$length_bp, evidence$plot_taxon, sum), decreasing = TRUE)
taxa <- head(setdiff(names(taxon_span), names(special)), 10)
shown <- c(taxa, names(special))
evidence$display_taxon <- ifelse(evidence$plot_taxon %in% shown, evidence$plot_taxon, "other taxa")
palette <- c(setNames(grDevices::hcl.colors(length(taxa), "Dark 3"), taxa),
             special, "other taxa" = "#C98942")
has_depth <- any(!is.na(evidence$mean_depth))
evidence$y <- if (has_depth) evidence$mean_depth else evidence$length_bp
points <- evidence[is.finite(evidence$gc_pct) & is.finite(evidence$y), ]
visible_taxa <- intersect(names(palette), unique(points$display_taxon))
missing_gc <- sum(!is.finite(evidence$gc_pct))
missing_depth <- if (has_depth) sum(!is.finite(evidence$mean_depth)) else nrow(evidence)
zero_depth <- sum(evidence$mean_depth == 0, na.rm = TRUE)
y_title <- if (has_depth) "Mean read depth" else "Sequence length (bp)"
transform <- if (has_depth) scales::pseudo_log_trans(base = 10, sigma = .1) else scales::log10_trans()
positive_depth <- evidence$mean_depth[is.finite(evidence$mean_depth) & evidence$mean_depth > 0]
depth_breaks <- if (has_depth && length(positive_depth) && max(positive_depth) >= 1 &&
                    max(positive_depth) / min(positive_depth) >= 10) {
    c(0, 10^seq(max(-4, floor(log10(min(positive_depth)))), ceiling(log10(max(positive_depth)))))
} else {
    waiver()
}
y_range <- range(evidence$y[is.finite(evidence$y)])
if (y_range[1] == y_range[2]) {
    y_range[2] <- y_range[2] + max(1, .1 * y_range[2])
}
display_range <- transform$transform(y_range)
display_range <- transform$inverse(display_range + c(-1, 1) * .05 * diff(display_range))
if (has_depth) display_range[1] <- max(0, display_range[1])
depth_scale <- function() scale_y_continuous(trans = transform, breaks = depth_breaks,
                                             labels = scales::label_number())
gc_scale <- function() scale_x_continuous(breaks = seq(0, 100, 20))
format_span <- scales::label_number(scale_cut = scales::cut_short_scale())
base <- theme_bw(base_size = 11) + theme(legend.position = "bottom")
scatter <- ggplot(points, aes(gc_pct, y, color = display_taxon)) +
    geom_point(alpha = .55, size = .45) +
    geom_point(aes(size = length_bp), alpha = .62) +
    gc_scale() + depth_scale() + scale_color_manual(values = palette, limits = visible_taxa) +
    coord_cartesian(xlim = c(0, 100), ylim = display_range) +
    scale_size_area(max_size = 9, labels = format_span,
                    breaks = unique(signif(max(evidence$length_bp) * c(.25, .5, 1), 2))) +
    guides(color = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1))) +
    labs(x = "GC (%)", y = y_title, color = "Taxonomic group", size = "Sequence length (bp)") +
    base
if (!nrow(points)) {
    scatter <- scatter + annotate("text", x = 50, y = mean(y_range),
                                  label = "No sequences have measurable GC and a plottable value")
}
gc_margin <- ggplot(points, aes(gc_pct, weight = length_bp, fill = display_taxon)) +
    geom_histogram(binwidth = 2, boundary = 0, alpha = .75, show.legend = FALSE) +
    gc_scale() + scale_y_continuous(labels = format_span) +
    coord_cartesian(xlim = c(0, 100)) +
    scale_fill_manual(values = palette, limits = visible_taxa) + labs(x = NULL, y = "Span (bp)") +
    base + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
depth_margin <- ggplot(points, aes(y = y, weight = length_bp, fill = display_taxon)) +
    geom_histogram(orientation = "y", bins = 30, alpha = .75, show.legend = FALSE) +
    depth_scale() + scale_x_continuous(labels = format_span) +
    coord_cartesian(ylim = display_range) +
    scale_fill_manual(values = palette, limits = visible_taxa) + labs(x = "Span (bp)", y = NULL) +
    base + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
blob <- gc_margin + plot_spacer() + scatter + depth_margin + guide_area() +
    plot_layout(design = "AB\nCD\nEE", widths = c(4, 1.35),
                heights = c(1.2, 4, 1.5), guides = "collect") +
    plot_annotation(title = if (has_depth) "GC and read coverage by taxonomic group" else "GC and sequence length",
                    caption = sprintf("Plotted %s of %s contigs  |  No measurable GC: %s  |  Zero depth: %s  |  Depth unavailable: %s",
                                      nrow(points), nrow(evidence), missing_gc, zero_depth, missing_depth))
ggsave(file.path(args[2], "blob.png"), blob, width = 10.5, height = 8.2, dpi = 180, bg = "white")

ordered <- evidence[order(evidence$display_taxon, -evidence$length_bp), ]
ordered$evidence <- ordered$display_taxon
ordered$rank <- ave(ordered$length_bp, ordered$evidence, FUN = seq_along)
ordered$cumulative_bp <- ave(ordered$length_bp, ordered$evidence, FUN = cumsum)
curves <- rbind(ordered[c("rank", "cumulative_bp", "evidence")],
                data.frame(rank = 0, cumulative_bp = 0, evidence = unique(ordered$evidence)))
use_mb <- max(curves$cumulative_bp) >= 1e6
divisor <- if (use_mb) 1e6 else 1e3
rank_breaks <- if (max(curves$rank) <= 10) seq(0, max(curves$rank)) else waiver()
cumulative <- ggplot(curves, aes(rank, cumulative_bp / divisor, color = evidence)) +
    geom_step() + geom_point(size = .7) +
    scale_x_continuous(breaks = rank_breaks) +
    scale_color_manual(values = palette) +
    labs(x = "Sequences ranked by length within each evidence group",
         y = if (use_mb) "Cumulative sequence (Mb)" else "Cumulative sequence (kb)",
         color = "Taxonomic group") +
    theme_bw(base_size = 12)
ggsave(file.path(args[2], "cumulative.png"), cumulative, width = 8, height = 5, dpi = 180, bg = "white")
