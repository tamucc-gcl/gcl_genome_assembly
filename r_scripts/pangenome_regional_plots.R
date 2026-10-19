# Presentation of audited regional tables; no graph computation.
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop('Expected sharing table, group context, chromosome presence and output directory')
out <- args[[4]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tiers <- c('core', 'softcore', 'shell', 'cloud')
palette <- c(core = '#0072B2', softcore = '#D55E00', shell = '#009E73', cloud = '#CC79A7')
data <- read_tsv(args[[1]], show_col_types = FALSE) |> filter(class %in% tiers) |>
  mutate(haplotype = replace_na(haplotype, ''), chromosome = replace_na(chromosome, ''),
         class = factor(class, levels = tiers))
context <- read_tsv(args[[2]], show_col_types = FALSE)
presence <- read_tsv(args[[3]], show_col_types = FALSE)
if (anyDuplicated(context[c('haplotype', 'chromosome')]) ||
    anyDuplicated(presence[c('haplotype', 'chromosome')])) stop('Duplicate regional context')
scoped <- data |> filter(scope == 'chromosome_haplotype') |>
  left_join(context, by = c('haplotype', 'chromosome'))
if (any(is.na(scoped$kind))) stop('Missing final-assembly context for a regional group')
# Expected chromosome-scale composites missing from CLIP stay unavailable, never zero.
missing_groups <- context |> filter(graph_group_present == 'no') |>
  crossing(class = factor(tiers, levels = tiers)) |>
  mutate(scope = 'chromosome_haplotype', bp = NA_real_, unit_distinct_bp = NA_real_, percent = NA_real_)
scoped <- bind_rows(scoped, missing_groups)
large_composites <- scoped |> filter(chromosome_scale_composite == 'yes') |> pull(chromosome) |> unique()
other_composites <- scoped |> filter(kind == 'composite', chromosome_scale_composite == 'no')
other_composites |> write_tsv(file.path(out, 'other_composite_sharing.tsv'))
presence |> write_tsv(file.path(out, 'chromosome_presence.tsv'))
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
    ggplot(rows, aes(.data[[variable]], label, fill = class)) + geom_col(width = 0.8, na.rm = TRUE) +
      geom_text(data = rows |> filter(is.na(percent)) |> distinct(label),
                aes(x = 0, y = label, label = 'No matching CLIP path group'),
                inherit.aes = FALSE, hjust = 0, size = 3) +
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
has_composites <- composition(scoped |> filter(chromosome_scale_composite == 'yes') |>
            mutate(composite_haplotype = str_c(haplotype, ' / ', str_remove(chromosome, '^composite:'))),
            'composite_haplotype', 'Chromosome-scale composite scaffolds by haplotype', 'composite_composition')

heatmaps <- function(rows, name, title, mark_composites = FALSE) {
  if (!nrow(rows)) return(FALSE)
  haplotypes <- data |> filter(scope == 'haplotype') |> pull(haplotype) |> unique() |> str_sort(numeric = TRUE)
  all_names <- str_sort(unique(rows$chromosome), numeric = TRUE)
  chromosomes <- c(all_names[!str_starts(all_names, 'composite:')], all_names[str_starts(all_names, 'composite:')])
  # Explicit NA cells: absent combinations must not become biological zeros.
  grid <- expand_grid(haplotype = haplotypes, chromosome = chromosomes, class = factor(tiers, levels = tiers)) |>
    left_join(rows |> select(haplotype, chromosome, class, percent), by = c('haplotype', 'chromosome', 'class')) |>
    left_join(presence |> select(haplotype, chromosome, status), by = c('haplotype', 'chromosome')) |>
    left_join(context |> select(haplotype, chromosome, graph_group_present), by = c('haplotype', 'chromosome')) |>
    mutate(marker = if (mark_composites) case_when(status == 'composite_only' ~ 'C',
               status == 'separate_and_composite' ~ '+', status == 'other_composite_only' ~ 'f',
               status == 'separate_and_composite_missing' ~ '+!', status == 'composite_missing_from_graph' ~ '!',
               graph_group_present == 'no' ~ '!', TRUE ~ '') else if_else(graph_group_present == 'no', '!', '', missing = ''),
           haplotype = factor(haplotype, levels = haplotypes), chromosome = factor(chromosome, levels = rev(chromosomes)))
  panels <- map(tiers, function(tier) {
    ggplot(grid |> filter(class == tier), aes(haplotype, chromosome, fill = percent)) +
      geom_tile(color = 'white', linewidth = 0.2) +
      geom_text(aes(label = marker), size = 3.1, color = 'black') +
      scale_fill_gradient(low = 'white', high = palette[[tier]], limits = c(0, 100), na.value = '#BBBBBB') +
      scale_y_discrete(labels = function(x) str_replace(x, '^composite:', 'Composite: ')) +
      labs(title = str_to_title(tier), x = NULL, y = NULL, fill = '% of unit') +
      theme(axis.text.x = element_text(angle = 55, hjust = 1), panel.grid = element_blank())
  })
  figure <- wrap_plots(panels, ncol = 2) + plot_annotation(title = title,
    subtitle = if (mark_composites) 'Whole-scaffold naming groups; composite sequence is counted in its own row' else
      'Global cohort sharing thresholds; grey = no reported group, not biological absence',
    caption = if (mark_composites) 'C = in chromosome-scale composite; + = separate pieces plus composite; ! = expected composite has no matching CLIP path group\nf = only smaller/fragmented composite labels. Grey without a symbol = no reported group. Symbols do not assign a percentage.' else
      '! = present in final assembly, no matching CLIP path group; no percentage is available')
  save_plot(name, figure, width = max(13, length(haplotypes) * 0.7),
            height = max(9, length(chromosomes) * 0.5 + 4))
  TRUE
}
main_rows <- scoped |> filter(kind != 'composite' | chromosome_scale_composite == 'yes')
heatmaps(main_rows, 'chromosome_haplotype_sharing', 'Sharing by chromosome naming group and haplotype', TRUE)
heatmaps(scoped |> filter(chromosome_scale_composite == 'yes'),
         'composite_haplotype_sharing', 'Chromosome-scale composite scaffold groups')
write_lines(c('# Regional sharing figures', '',
  '<details><summary>Haplotype and chromosome composition</summary>', '',
  '![Haplotype composition](haplotype_composition.png)', '',
  '![Chromosome composition](chromosome_composition.png)', '',
  '![Chromosome-by-haplotype sharing](chromosome_haplotype_sharing.png)', '',
  'Rows are whole-scaffold naming groups, not a chromosome-by-coordinate painting. A chr4 row contains separately named chr4 pieces; a chr4+chr14 scaffold is counted only in its composite row. Consequently chr4/chr7 percentages may describe only their separate pieces, and grey chr12/chr14 cells do not mean those sequences are absent.', '',
  '`C` marks a chromosome named in a chromosome-scale composite; `+` marks separate pieces plus a chromosome-scale composite. `f` means only smaller/fragmented composite labels were reported. These symbols describe representation, not measured sharing percentages for the missing row.', '',
  '`!` flags a chromosome-scale composite present in the final assembly but absent under its expected name in CLIP regional groups (`+!` also has separately named pieces). This is missing graph representation requiring clipping/path checks, not biological absence. Missing measurements remain blank.', '',
  'Coverage thresholds use the entire cohort. Each unit counts distinct graph sequence once. Chromosome unions can overlap and must not be summed. Unplaced sequence is retained; private overlaps the tier classification and is not added as a fifth tier.', '',
  '[Chromosome representation table](chromosome_presence.tsv)', '',
  '[Source table](../regional_sharing.tsv) | [Audit](../regional_audit.json)', '', '</details>', '',
  if (has_composites) c('<details><summary>Chromosome-scale composite scaffold groups</summary>', '',
    'These groups contain scaffolds selected as chromosome-scale in a chromosome-scale final assembly. Composite labels are not split, double-counted across their named chromosomes, or interpreted as confirmed fusions.', '',
    '![Composite composition](composite_composition.png)', '',
    '![Composite-by-haplotype sharing](composite_haplotype_sharing.png)', '', '</details>') else '',
  '<details><summary>Other composite labels (supplementary)</summary>', '',
  str_glue('{n_distinct(other_composites$chromosome)} other composite naming groups are retained in the tables. Small scaffolds and assemblies without chromosome-scale structure are excluded from the main composite figures; names alone are not evidence of chromosome fusions.'), '',
  '[All other composite measurements](other_composite_sharing.tsv). No sequence is removed from the graph or global sharing counts.', '', '</details>'),
  file.path(out, 'regional_figures.md'))
packages <- c('tidyverse', 'ggplot2', 'patchwork')
tibble(tool = c('R', packages), version = c(as.character(getRversion()), map_chr(packages, ~ as.character(packageVersion(.x))))) |>
  mutate(process = 'PANGENOME_REGIONAL_PLOTS', .before = 1) |> write_tsv(file.path(out, 'plot_versions.tsv'))
