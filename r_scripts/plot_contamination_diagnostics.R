library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
evidence <- read.delim(args[1], check.names = FALSE)
dir.create(args[2], showWarnings = FALSE)
palette <- c("host-compatible" = "#0072B2", "non-host" = "#D55E00",
             "mixed" = "#AA4499", "unresolved" = "#777777")

has_depth <- any(!is.na(evidence$mean_depth))
if (has_depth) {
    blob <- ggplot(evidence, aes(100 * gc_fraction, mean_depth,
                                  color = evidence, size = length_bp))
    blob <- blob + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
        labs(y = "Mean read depth", title = "GC and read coverage")
} else {
    blob <- ggplot(evidence, aes(100 * gc_fraction, length_bp, color = evidence))
    blob <- blob + scale_y_log10() +
        labs(y = "Sequence length (bp)", title = "GC and sequence length",
             subtitle = "Read coverage was unavailable")
}
blob <- blob + geom_point(alpha = .75, na.rm = TRUE) +
    scale_color_manual(values = palette) +
    scale_size_area(max_size = 9, labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = "GC (%)", color = "Taxonomic evidence", size = "Sequence length (bp)",
         caption = sprintf("Sequences without measurable GC: %s", sum(is.na(evidence$gc_fraction)))) +
    theme_bw(base_size = 12)
ggsave(file.path(args[2], "blob.png"), blob, width = 8, height = 5, dpi = 180, bg = "white")

ordered <- evidence[order(evidence$evidence, -evidence$length_bp), ]
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
         color = "Taxonomic evidence") +
    theme_bw(base_size = 12)
ggsave(file.path(args[2], "cumulative.png"), cumulative, width = 8, height = 5, dpi = 180, bg = "white")
