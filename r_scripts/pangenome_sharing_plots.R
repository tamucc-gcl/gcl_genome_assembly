# Plot validated Panacus outputs; presentation never recomputes growth.
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Expected input and output directories')
input <- args[[1]]
output <- args[[2]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
units <- c('haplotype', 'individual')
tiers <- c('core', 'softcore', 'shell', 'cloud')
palette <- c(core = '#0072B2', softcore = '#D55E00', shell = '#009E73', cloud = '#CC79A7')
theme_set(theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank()))
read_panacus <- function(path) {
  lines <- read_lines(path) |> discard(~ str_starts(.x, '#') || !nzchar(str_trim(.x)))
  fields <- str_split(lines, fixed('\t'))
  numeric <- map_lgl(fields, ~ str_detect(.x[[1]], '^[0-9]+$'))
  headers <- fields[!numeric] |> set_names(map_chr(fields[!numeric], 1)) |> map(~ .x[-1])
  if (!any(numeric)) stop('No numeric rows: ', path)
  data <- read_tsv(I(str_c(lines[numeric], collapse = '\n')), col_names = FALSE,
                   col_types = cols(.default = col_double()), show_col_types = FALSE)
  if (nrow(problems(data)) || any(!is.finite(as.matrix(data)))) stop('Invalid table: ', path)
  list(headers = headers, data = data)
}
save_panels <- function(name, panels) {
  figure <- wrap_plots(panels, ncol = 2, guides = 'collect') & theme(legend.position = 'bottom')
  walk(c('png', 'pdf'), function(ext) {
    ggsave(file.path(output, str_c(name, '.', ext)), figure, width = 13, height = 6.5,
           dpi = 180, bg = 'white')
  })
}
sharing <- read_tsv(file.path(input, 'sharing_summary.tsv'), show_col_types = FALSE)
thresholds <- sharing |> distinct(unit, denominator, core_min, softcore_min, shell_min)
if (anyDuplicated(thresholds$unit) || !setequal(thresholds$unit, units)) stop('Inconsistent thresholds')
threshold_label <- function(unit_name) {
  row <- thresholds |> filter(unit == unit_name)
  str_glue('N={row$denominator}; core >={row$core_min}, softcore >={row$softcore_min}, shell >={row$shell_min}')
}
save_panels('growth', map(units, function(unit_name) {
  table <- read_panacus(file.path(input, str_c(unit_name, '.growth.tsv')))
  cov <- as.numeric(table$headers$coverage)
  quo <- as.numeric(table$headers$quorum)
  if (length(cov) != ncol(table$data) - 1L || length(quo) != length(cov)) stop('Growth header mismatch')
  labels <- str_glue('Coverage >={cov}; quorum {100 * quo}%')
  keys <- str_c('curve_', seq_along(labels))
  data <- table$data |> set_names(c('cohort_size', keys)) |>
    pivot_longer(-cohort_size, names_to = 'curve', values_to = 'bp') |>
    mutate(curve = factor(curve, levels = keys, labels = labels))
  ggplot(data, aes(cohort_size, bp / 1e6, color = curve, linetype = curve)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    scale_x_continuous(breaks = sort(unique(data$cohort_size))) +
    scale_y_continuous(limits = c(0, NA), labels = scales::label_comma()) +
    labs(title = str_c(str_to_title(unit_name), ' accumulation'),
         x = str_c('Number of ', unit_name, 's'), y = 'Distinct graph sequence (Mb)',
         color = 'Growth criterion', linetype = 'Growth criterion')
}))
save_panels('sharing_histogram', map(units, function(unit_name) {
  data <- read_panacus(file.path(input, str_c(unit_name, '.hist.tsv')))$data
  if (ncol(data) != 2L) stop('Expected one histogram column')
  data <- data |> set_names(c('coverage', 'bp')) |> filter(coverage > 0) |>
    mutate(unit = unit_name) |> left_join(thresholds, by = 'unit') |>
    mutate(class = case_when(coverage >= core_min ~ 'core', coverage >= softcore_min ~ 'softcore',
                            coverage >= shell_min ~ 'shell', TRUE ~ 'cloud'),
           class = factor(class, levels = tiers))
  observed <- data |> group_by(class, .drop = FALSE) |> summarise(bp = sum(bp), .groups = 'drop')
  expected <- sharing |> filter(unit == unit_name, class %in% tiers) |>
    transmute(class = factor(class, levels = tiers), expected_bp = bp)
  check <- full_join(observed, expected, by = 'class')
  if (nrow(check) != 4L || anyNA(check) || any(check$bp != check$expected_bp)) stop('Histogram tiers disagree with partition')
  ggplot(data, aes(coverage, bp / 1e6, fill = class)) + geom_col(width = 0.8) +
    scale_fill_manual(values = palette, limits = tiers, drop = FALSE) +
    scale_x_continuous(breaks = data$coverage) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.05))) +
    labs(title = str_c(str_to_title(unit_name), ' sharing'), subtitle = threshold_label(unit_name),
         x = str_c('Number of ', unit_name, 's containing sequence'), y = 'Distinct graph sequence (Mb)', fill = 'Sharing tier')
}))
save_panels('sharing_partition', map(units, function(unit_name) {
  data <- sharing |> filter(unit == unit_name, class %in% tiers) |> mutate(class = factor(class, levels = tiers))
  ggplot(data, aes(class, percent_represented_graph_bp, fill = class)) + geom_col(width = 0.8) +
    scale_fill_manual(values = palette, limits = tiers, drop = FALSE) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
    labs(title = str_c(str_to_title(unit_name), ' partition'), subtitle = threshold_label(unit_name),
         x = NULL, y = '% represented graph bp', fill = 'Sharing tier')
}))
write_lines(c('# CLIP sharing figures', '', '<details><summary>Growth and sharing</summary>', '',
  '![Growth curves](growth.png)', '',
  'Panacus values are plotted directly, without openness fits or invented uncertainty bands. Growth quorum is evaluated at each cohort size and differs from the all-but-one softcore partition threshold.', '',
  '![Sharing histogram](sharing_histogram.png)', '', '![Sharing partition](sharing_partition.png)', '',
  'Histogram and partition colors use the same thresholds and precedence: core, softcore, shell, cloud. Private sequence is not an additional disjoint tier. Values describe graph representation, not biological novelty.', '',
  'Source tables: [sharing summary](../sharing/sharing_summary.tsv), [haplotype growth](../sharing/haplotype.growth.tsv), [individual growth](../sharing/individual.growth.tsv).', '', '</details>'),
  file.path(output, 'sharing_figures.md'))
packages <- c('tidyverse', 'ggplot2', 'patchwork')
tibble(tool = c('R', packages), version = c(as.character(getRversion()), map_chr(packages, ~ as.character(packageVersion(.x))))) |>
  mutate(process = 'PANGENOME_SHARING_PLOTS', .before = 1) |> write_tsv(file.path(output, 'versions.tsv'))
