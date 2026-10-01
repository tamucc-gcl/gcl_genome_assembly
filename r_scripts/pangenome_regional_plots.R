# Presentation of audited regional tables; no graph computation.
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Expected regional_sharing.tsv and output directory')
out <- args[[2]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tiers <- c('core', 'softcore', 'shell', 'cloud')
palette <- c(core = '#0072B2', softcore = '#D55E00', shell = '#009E73', cloud = '#CC79A7')
data <- read_tsv(args[[1]], show_col_types = FALSE) |> filter(class %in% tiers) |>
  mutate(haplotype = replace_na(haplotype, ''), chromosome = replace_na(chromosome, ''),
         class = factor(class, levels = tiers))
checks <- data |> group_by(scope, haplotype, chromosome) |>
  summarise(rows = n(), classes = n_distinct(class), bp_sum = sum(bp),
            total = first(unit_distinct_bp), totals = n_distinct(unit_distinct_bp), .groups = 'drop')
if (any(checks$rows != 4 | checks$classes != 4 | checks$totals != 1 | checks$bp_sum != checks$total) ||
    any(!is.finite(data$bp) | data$bp < 0) ||
    any(abs(data$percent[data$unit_distinct_bp > 0] -
            100 * data$bp[data$unit_distinct_bp > 0] / data$unit_distinct_bp[data$unit_distinct_bp > 0]) > 1e-5)) {
  stop('Invalid regional partition table')
}
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank()))
save_plot <- function(name, figure, width = 13, height = 7) {
  walk(c('png', 'pdf'), ~ ggsave(file.path(out, str_c(name, '.', .x)), figure,
                               width = width, height = height, dpi = 180, bg = 'white'))
}
composition <- function(rows, key, title, name) {
  if (!nrow(rows)) return(FALSE)
  levels <- str_sort(unique(rows[[key]]), numeric = TRUE)
  rows <- rows |> mutate(label = factor(.data[[key]], levels = rev(levels)))
  panel <- function(variable, axis) {
    ggplot(rows, aes(.data[[variable]], label, fill = class)) + geom_col(width = 0.8) +
      scale_fill_manual(values = palette, limits = tiers, drop = FALSE) +
      scale_x_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.03))) +
      labs(x = axis, y = NULL, fill = 'Sharing tier')
  }
  rows <- rows |> mutate(mb = bp / 1e6)
  figure <- (panel('mb', 'Distinct graph sequence (Mb)') |
             panel('percent', '% of this reporting unit')) +
    plot_layout(guides = 'collect') + plot_annotation(title = title)
  figure <- figure & theme(legend.position = 'bottom')
  save_plot(name, figure, height = max(5, length(levels) * 0.28 + 2))
  TRUE
}
composition(data |> filter(scope == 'haplotype'), 'haplotype',
            'Haplotype sharing composition', 'haplotype_composition')
composition(data |> filter(scope == 'chromosome', !str_starts(chromosome, 'composite:')),
            'chromosome', 'Chromosome unions and unplaced sequence (not additive across rows)', 'chromosome_composition')
has_composites <- composition(data |> filter(scope == 'chromosome', str_starts(chromosome, 'composite:')),
            'chromosome', 'Composite scaffold groups (kept separate from chromosomes)', 'composite_composition')

heatmaps <- function(rows, name, title) {
  if (!nrow(rows)) return(FALSE)
  haplotypes <- data |> filter(scope == 'haplotype') |> pull(haplotype) |> unique() |> str_sort(numeric = TRUE)
  chromosomes <- str_sort(unique(rows$chromosome), numeric = TRUE)
  # Explicit NA cells: absent combinations must not become biological zeros.
  grid <- expand_grid(haplotype = haplotypes, chromosome = chromosomes, class = factor(tiers, levels = tiers)) |>
    left_join(rows |> select(haplotype, chromosome, class, percent), by = c('haplotype', 'chromosome', 'class')) |>
    mutate(haplotype = factor(haplotype, levels = haplotypes), chromosome = factor(chromosome, levels = rev(chromosomes)))
  panels <- map(tiers, function(tier) {
    ggplot(grid |> filter(class == tier), aes(haplotype, chromosome, fill = percent)) +
      geom_tile(color = 'white', linewidth = 0.2) +
      scale_fill_gradient(low = 'white', high = palette[[tier]], limits = c(0, 100), na.value = '#BBBBBB') +
      labs(title = str_to_title(tier), x = NULL, y = NULL, fill = '% of unit') +
      theme(axis.text.x = element_text(angle = 55, hjust = 1), panel.grid = element_blank())
  })
  figure <- wrap_plots(panels, ncol = 2) + plot_annotation(title = title,
    subtitle = 'Global cohort sharing thresholds; grey = no reported combination, not zero')
  save_plot(name, figure, width = max(13, length(haplotypes) * 0.7),
            height = max(9, length(chromosomes) * 0.5 + 4))
  TRUE
}
heatmaps(data |> filter(scope == 'chromosome_haplotype', !str_starts(chromosome, 'composite:')),
         'chromosome_haplotype_sharing', 'Sharing composition by chromosome and haplotype')
heatmaps(data |> filter(scope == 'chromosome_haplotype', str_starts(chromosome, 'composite:')),
         'composite_haplotype_sharing', 'Sharing composition of composite scaffold groups')
write_lines(c('# Regional sharing figures', '',
  '<details><summary>Haplotype and chromosome composition</summary>', '',
  '![Haplotype composition](haplotype_composition.png)', '',
  '![Chromosome composition](chromosome_composition.png)', '',
  '![Chromosome-by-haplotype sharing](chromosome_haplotype_sharing.png)', '',
  'Coverage thresholds use the entire cohort. Each unit counts distinct graph sequence once. Chromosome unions can overlap and must not be summed. Grey cells are unavailable combinations. Unplaced sequence is retained; private overlaps the tier classification and is not added as a fifth tier.', '',
  '[Source table](../regional_sharing.tsv) · [Audit](../regional_audit.json)', '', '</details>', '',
  if (has_composites) c('<details><summary>Composite scaffold groups</summary>', '',
    'Composite labels are not split or reassigned to individual chromosomes.', '',
    '![Composite composition](composite_composition.png)', '',
    '![Composite-by-haplotype sharing](composite_haplotype_sharing.png)', '', '</details>') else ''),
  file.path(out, 'regional_figures.md'))
packages <- c('tidyverse', 'ggplot2', 'patchwork')
tibble(tool = c('R', packages), version = c(as.character(getRversion()), map_chr(packages, ~ as.character(packageVersion(.x))))) |>
  mutate(process = 'PANGENOME_REGIONAL_PLOTS', .before = 1) |> write_tsv(file.path(out, 'plot_versions.tsv'))
