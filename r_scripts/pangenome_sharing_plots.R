# Presentation only: plot validated Panacus tables, never recompute growth.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: pangenome_sharing_plots.R input_directory output_directory')
input <- args[1]; output <- args[2]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
read_panacus <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl('^#', lines) & nzchar(trimws(lines))]
  fields <- strsplit(lines, '\t', fixed = TRUE)
  numeric <- vapply(fields, function(x) grepl('^[0-9]+$', x[1]), logical(1))
  headers <- fields[!numeric]
  names(headers) <- vapply(headers, `[`, character(1), 1)
  data <- do.call(rbind, lapply(fields[numeric], as.numeric))
  if (is.null(data) || any(!is.finite(data))) stop('Invalid table: ', path)
  list(headers = headers, data = data)
}
draw <- function(name, fn) {
  png(file.path(output, paste0(name, '.png')), width = 1800, height = 900, res = 150, type = 'cairo')
  tryCatch(fn(), finally = dev.off())
  pdf(file.path(output, paste0(name, '.pdf')), width = 12, height = 6, useDingbats = FALSE)
  tryCatch(fn(), finally = dev.off())
}
units <- c('haplotype', 'individual')
colors <- c('#0072B2', '#D55E00', '#009E73', '#CC79A7')
growth <- lapply(units, function(unit) read_panacus(file.path(input, paste0(unit, '.growth.tsv'))))
hist <- lapply(units, function(unit) read_panacus(file.path(input, paste0(unit, '.hist.tsv'))))
summary <- read.delim(file.path(input, 'sharing_summary.tsv'), check.names = FALSE)
draw('growth', function() {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  for (i in seq_along(units)) {
    table <- growth[[i]]; d <- table$data
    cov <- as.numeric(table$headers$coverage[-1]); quo <- as.numeric(table$headers$quorum[-1])
    if (length(cov) != ncol(d) - 1L || length(quo) != length(cov)) stop('Growth header mismatch')
    matplot(d[, 1], d[, -1, drop = FALSE] / 1e6, type = 'o', pch = seq_along(cov),
            col = rep(colors, length.out = length(cov)), lty = seq_along(cov),
            xlab = paste('Number of', paste0(units[i], 's')), ylab = 'Distinct graph sequence (Mb)',
            main = paste(tools::toTitleCase(units[i]), 'accumulation'), ylim = c(0, max(d[, -1]) / 1e6))
    labels <- paste0('Coverage >= ', cov, '; quorum ', 100 * quo, '%')
    legend('topleft', legend = labels, col = rep(colors, length.out = length(cov)),
           lty = seq_along(cov), pch = seq_along(cov), bty = 'n', cex = 0.75)
  }
})
draw('sharing_histogram', function() {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  for (i in seq_along(units)) {
    d <- hist[[i]]$data; d <- d[d[, 1] > 0, , drop = FALSE]
    barplot(d[, 2] / 1e6, names.arg = d[, 1], col = colors[1], border = NA,
            xlab = paste('Number of', paste0(units[i], 's'), 'containing sequence'),
            ylab = 'Distinct graph sequence (Mb)', main = paste(tools::toTitleCase(units[i]), 'sharing'))
  }
})
tiers <- c('core', 'softcore', 'shell', 'cloud')
draw('sharing_partition', function() {
  par(mfrow = c(1, 2), mar = c(6, 5, 4, 1))
  for (unit in units) {
    rows <- summary[summary$unit == unit & summary$class %in% tiers, ]
    rows <- rows[match(tiers, rows$class), ]
    if (anyNA(rows$bp)) stop('Missing sharing tier')
    barplot(rows$percent_represented_graph_bp, names.arg = tiers, col = colors, border = NA,
            ylab = '% represented graph bp', ylim = c(0, 100), main = paste(tools::toTitleCase(unit), 'partition'))
    mtext(sprintf('N=%s; minima: core %s, softcore %s, shell %s', rows$denominator[1],
                  rows$core_min[1], rows$softcore_min[1], rows$shell_min[1]), side = 1, line = 4, cex = 0.75)
  }
})
report <- c('# CLIP sharing figures', '',
  'Panacus growth values are plotted directly. No fitted openness classification or uncertainty band is added.', '',
  '<details><summary>Growth and sharing</summary>', '',
  '![Growth curves](growth.png)', '',
  'Growth quorum is evaluated at each cohort size. A 90% quorum requires 5/5 individuals in a five-individual cohort; this differs from the summary tier defined as all-but-one.', '',
  '![Sharing histogram](sharing_histogram.png)', '',
  '![Sharing partition](sharing_partition.png)', '',
  'Partition tiers are disjoint and use the displayed minimum counts. Private sequence is not added as another tier. Values describe graph sequence representation, not gene presence/absence or biological novelty.', '',
  'Source tables: [sharing summary](../sharing/sharing_summary.tsv), [haplotype growth](../sharing/haplotype.growth.tsv), [individual growth](../sharing/individual.growth.tsv).', '',
  '</details>', '')
writeLines(report, file.path(output, 'sharing_figures.md'))
writeLines(c('process\ttool\tversion', paste('PANGENOME_SHARING_PLOTS', 'R', getRversion(), sep = '\t')),
           file.path(output, 'versions.tsv'))
