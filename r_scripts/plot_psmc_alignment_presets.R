#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]
d <- read.delim(file.path(root, "classification_summary.tsv"), check.names = FALSE)
d$preset <- factor(d$preset, levels = c("asm5", "asm10", "asm20"))
states <- c("no_PAF_coverage", "secondary_only", "incomplete_scoring_only", "below_10kb_only",
            "below_MAPQ5_only", "overlapping_segments_one_query", "overlapping_multiple_queries",
            "unique_10_to_50kb", "unique_long_alignment")
labels <- setNames(vapply(states, function(s) d$label[match(s, d$state)], character(1)), states)
d$state <- factor(d$state, levels = states)
colours <- setNames(c("#666666", "#AAAAAA", "#DDDDDD", "#CC6677", "#882255",
                     "#AA4499", "#332288", "#DDCC77", "#228833"), states)
p <- ggplot(d, aes(preset, percent, fill = state)) +
    geom_col(width = .7) + facet_wrap(~haplotype, nrow = 1) +
    scale_fill_manual(values = colours, labels = function(x) vapply(labels[x], function(s)
        paste(strwrap(s, 32), collapse = "\n"), character(1)), drop = FALSE) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, .02))) +
    labs(title = "Dipcall mask support with alternative assembly presets",
         subtitle = "Same reference scope and confidence rules; asm5 is the reused control",
         x = "Minimap2 preset", y = "% of selected nongap reference sequence", fill = NULL) +
    theme_bw(base_size = 12) + theme(panel.grid.major.x = element_blank(),
        legend.position = "bottom", plot.background = element_rect(fill = "white", colour = NA)) +
    guides(fill = guide_legend(ncol = 3, byrow = TRUE))
ggsave(file.path(root, "preset_comparison.png"), p, width = 11, height = 8, dpi = 160, bg = "white")
ggsave(file.path(root, "preset_comparison.pdf"), p, width = 11, height = 8, bg = "white")
write.table(data.frame(tool = c("R", "ggplot2"), version = c(as.character(getRversion()),
    as.character(packageVersion("ggplot2")))), file.path(root, "plot_versions.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
