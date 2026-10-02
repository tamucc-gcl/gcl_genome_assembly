suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
source_dir <- args[[1]]
out <- args[[2]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
coords <- read_tsv(file.path(source_dir, 'out.filtered.coords'),
  col_names = c('rstart', 'rend', 'qstart', 'qend', 'rlen', 'qlen', 'identity',
                'rdir', 'qdir', 'ref', 'query'), show_col_types = FALSE)
calls <- read_tsv(file.path(source_dir, 'syri.out'), col_types = cols(.default = col_character()),
  col_names = c('ref', 'rstart', 'rend', 'refseq', 'altseq', 'query', 'qstart', 'qend',
                'id', 'parent', 'type', 'copy_status'), na = '-', show_col_types = FALSE) |>
  mutate(across(c(rstart, rend, qstart, qend), as.numeric))
lengths <- read_tsv(file.path(out, 'sequence_lengths.tsv'), show_col_types = FALSE)
stopifnot(n_distinct(coords$ref) == 1, n_distinct(coords$query) == 1)
types <- c('INV', 'TRANS', 'INVTR', 'DUP', 'INVDP', 'CPG', 'CPL', 'TDM')
events <- calls |> filter(type %in% types) |>
  mutate(reference_span = abs(rend - rstart) + 1,
         query_span = abs(qend - qstart) + 1,
         rank_span = pmax(reference_span, query_span, na.rm = TRUE))
write_tsv(events |> select(-refseq, -altseq), file.path(out, 'candidate_events.tsv'))
counts <- events |> group_by(type) |>
  summarise(regions = n(), median_span = median(rank_span), max_span = max(rank_span), .groups = 'drop')
write_tsv(counts, file.path(out, 'event_counts.tsv'))
# Union intervals before measuring coverage: repeated matches must not inflate it.
union_intervals <- function(start, end) {
  low <- pmin(start, end) - 1
  high <- pmax(start, end)
  x <- tibble(start = low, end = high) |> arrange(start, end)
  if (nrow(x) == 0) return(x)
  x |> mutate(group = cumsum(start > lag(cummax(end), default = -Inf))) |>
    group_by(group) |> summarise(start = min(start), end = max(end), .groups = 'drop') |>
    select(-group)
}
coverage <- map_dfr(c('reference', 'query'), function(role) {
  intervals <- if (role == 'reference') union_intervals(coords$rstart, coords$rend) else
    union_intervals(coords$qstart, coords$qend)
  total <- lengths |> filter(.data$role == .env$role) |> pull(bases)
  stopifnot(length(total) == 1, max(intervals$end) <= total)
  edges <- seq(0, total, length.out = 101)
  tibble(role = role, start = head(edges, -1), end = tail(edges, -1)) |>
    mutate(aligned_bases = map2_dbl(start, end, function(a, b)
      sum(pmax(0, pmin(intervals$end, b) - pmax(intervals$start, a)))),
      fraction = aligned_bases / (end - start))
})
write_tsv(coverage, file.path(out, 'alignment_coverage_bins.tsv'))
coverage_totals <- coverage |> group_by(role) |>
  summarise(aligned_bases = sum(aligned_bases), bases = sum(end - start),
            aligned_fraction = aligned_bases / bases, .groups = 'drop')
write_tsv(coverage_totals, file.path(out, 'alignment_coverage.tsv'))
theme_set(theme_bw(base_size = 11))
dot <- ggplot(coords, aes(rstart / 1e6, qstart / 1e6,
                          xend = rend / 1e6, yend = qend / 1e6,
                          colour = if_else(qstart <= qend, 'Forward', 'Reverse'))) +
  geom_segment(linewidth = .7, alpha = .9) +
  scale_colour_manual(values = c(Forward = '#B34700', Reverse = '#0066A6')) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1))) +
  labs(x = 'Reference assembly (Mb)', y = 'Query assembly (Mb)', colour = 'Orientation',
       title = 'Filtered assembly-to-assembly alignment blocks')
sizes <- ggplot(events, aes(rank_span, fill = type)) +
  geom_histogram(bins = 35, show.legend = FALSE) + scale_x_log10() +
  facet_wrap(vars(type), scales = 'free_y') +
  labs(x = 'Larger reported reference/query span (bp; log scale)', y = 'Regions',
       title = 'Native SyRI classes; spans are not net sequence gain')
covplot <- ggplot(coverage, aes((start + end) / 2e6, fraction, colour = role)) +
  geom_line(linewidth = .8) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = 'Position in each assembly (Mb)', y = 'Union alignment coverage', colour = 'Assembly')
ggsave(file.path(out, 'overview.png'), ((dot + coord_equal()) / covplot) + plot_layout(heights = c(3, 1)),
       width = 10, height = 10, dpi = 180)
ggsave(file.path(out, 'event_sizes.png'), sizes, width = 12, height = 8, dpi = 180)
largest <- events |> group_by(type) |> slice_max(rank_span, n = 2, with_ties = FALSE) |> ungroup()
write_tsv(largest |> select(-refseq, -altseq), file.path(out, 'largest_candidates.tsv'))
walk(seq_len(nrow(largest)), function(i) {
  event <- largest[i, ]
  left <- max(0, min(event$rstart, event$rend) - 250000)
  right <- max(event$rstart, event$rend) + 250000
  bottom <- max(0, min(event$qstart, event$qend) - 250000)
  top <- max(event$qstart, event$qend) + 250000
  marked <- dot +
    geom_rect(data = event, inherit.aes = FALSE,
      aes(xmin = pmin(rstart, rend) / 1e6, xmax = pmax(rstart, rend) / 1e6,
          ymin = pmin(qstart, qend) / 1e6, ymax = pmax(qstart, qend) / 1e6),
      fill = NA, colour = 'black', linetype = 2, linewidth = .7)
  context <- marked + coord_fixed(xlim = c(left, right) / 1e6, ylim = c(bottom, top) / 1e6) +
    labs(title = 'Context: 250 kb flanks')
  pad <- max(10000, event$rank_span * .25)
  detail <- marked + coord_fixed(
    xlim = c(max(0, min(event$rstart, event$rend) - pad), max(event$rstart, event$rend) + pad) / 1e6,
    ylim = c(max(0, min(event$qstart, event$qend) - pad), max(event$qstart, event$qend) + pad) / 1e6) +
    labs(title = 'Candidate detail')
  plot <- (context + detail) + plot_layout(guides = 'collect') +
    plot_annotation(title = paste(event$type, event$id),
      subtitle = 'Assembly alignment blocks, not reads. Dashed box: reported reference/query intervals.')
  ggsave(file.path(out, sprintf('candidate_%02d.png', i)), plot, width = 14, height = 7, dpi = 180)
})
writeLines(c('# Chromosome-pair diagnostic report', '',
  'Exploratory alignment-based candidates, not validated biological events.',
  'Segments represent assembly-to-assembly matches, not individual reads. Blank space is not read depth.',
  'One chromosome pair cannot establish interchromosomal movement. No repeat annotation was supplied.',
  'Coverage is the union of filtered alignment spans, including internal gaps; it is not exact matching-base coverage.',
  'Bins use each assembly’s own coordinates. SyRI NOTAL statistics are a different measure.',
  'Alignment-child rows are excluded from event counts. Native CPG/CPL/TDM classes are retained.',
  'Reported spans can overlap and are not net gained/lost bases. Candidate zooms show two largest regions per class.',
  '', '![Alignment overview](overview.png)', '', '![Event sizes](event_sizes.png)', '',
  '## Candidate zooms', '',
  sprintf('![Candidate %02d](candidate_%02d.png)', seq_len(nrow(largest)), seq_len(nrow(largest)))),
  file.path(out, 'report.md'))
capture.output(sessionInfo(), file = file.path(out, 'R_session.txt'))
