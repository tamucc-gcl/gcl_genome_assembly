# Descriptive comparison of assembly haplotypes, not individual population inference.
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ape)
  library(ggrepel)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop('Expected similarity TSV, identity ledger and output directory')
out <- args[[3]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
identity <- read_tsv(args[[2]], show_col_types = FALSE, col_types = cols(.default = col_character()))
raw <- read_tsv(args[[1]], comment = '#', show_col_types = FALSE, name_repair = 'minimal')
if (names(raw)[[1]] != 'group') stop('Unexpected Panacus similarity header')
ids <- raw$group
if (anyDuplicated(ids) || anyDuplicated(identity$id) ||
    !setequal(ids, identity$id) || !setequal(names(raw)[-1], ids)) stop('Similarity identity mismatch')
matrix <- raw |> select(all_of(ids)) |> as.matrix()
storage.mode(matrix) <- 'double'
rownames(matrix) <- ids
if (nrow(matrix) != ncol(matrix) || any(!is.finite(matrix)) ||
    any(matrix < -1e-6 | matrix > 1 + 1e-6) ||
    max(abs(matrix - t(matrix))) > 1e-6 || max(abs(diag(matrix) - 1)) > 1e-6) {
  stop('Invalid similarity matrix')
}
# Canonical ledger ordering, not clustering order exported by Panacus.
ids <- sort(identity$id)
matrix <- matrix[ids, ids, drop = FALSE]
matrix <- (matrix + t(matrix)) / 2
matrix[matrix < 0] <- 0
matrix[matrix > 1] <- 1
diag(matrix) <- 1
distance <- 1 - matrix
n <- length(ids)
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank()))
save_plot <- function(name, plot, width = 10, height = 7) {
  walk(c('png', 'pdf'), ~ ggsave(file.path(out, str_c(name, '.', .x)), plot,
                               width = width, height = height, dpi = 180, bg = 'white'))
}
as_tibble(distance, rownames = 'id', .name_repair = 'minimal') |>
  write_tsv(file.path(out, 'haplotype_distance.tsv'))
heat <- as_tibble(matrix, rownames = 'id', .name_repair = 'minimal') |>
  pivot_longer(-id, names_to = 'other', values_to = 'similarity') |>
  mutate(id = factor(id, levels = rev(ids)), other = factor(other, levels = ids)) |>
  ggplot(aes(other, id, fill = similarity)) + geom_tile() +
  scale_fill_viridis_c(limits = c(0, 1)) + coord_equal() +
  labs(title = 'Haplotype sequence similarity', subtitle = 'Panacus bp-weighted Jaccard; CLIP graph',
       x = NULL, y = NULL, fill = 'Jaccard') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
save_plot('haplotype_similarity', heat, max(9, n * 0.45), max(7, n * 0.35))

# Classical PCoA via its defining double-centred squared-distance matrix.
# Do not label this PCA, silently correct distances, or discard negative eigenvalues from the audit.
center <- diag(n) - base::matrix(1 / n, n, n)
decomp <- eigen(-0.5 * center %*% (distance^2) %*% center, symmetric = TRUE)
tol <- max(1, max(abs(decomp$values))) * 1e-10
positive <- which(decomp$values > tol)
negative <- decomp$values[decomp$values < -tol]
positive_sum <- sum(decomp$values[positive])
negative_fraction <- if (sum(abs(decomp$values)) > 0) sum(abs(negative)) / sum(abs(decomp$values)) else 0
tibble(axis = seq_len(n), eigenvalue = decomp$values) |>
  write_tsv(file.path(out, 'pcoa_eigenvalues.tsv'))
coordinates <- base::matrix(0, n, 2)
for (axis in seq_len(min(2L, length(positive)))) {
  index <- positive[[axis]]
  coordinates[, axis] <- decomp$vectors[, index] * sqrt(decomp$values[[index]])
}
points <- tibble(id = ids, axis1 = coordinates[, 1], axis2 = coordinates[, 2]) |>
  left_join(identity |> select(id, sample, haplotype, is_reference), by = 'id')
write_tsv(points, file.path(out, 'pcoa_coordinates.tsv'))
axis_label <- function(axis) {
  if (length(positive) < axis) return(str_c('PCoA ', axis, ' (no positive axis)'))
  str_glue('PCoA {axis} ({round(100 * decomp$values[positive[[axis]]] / positive_sum, 1)}% positive inertia)')
}
ordination <- ggplot(points, aes(axis1, axis2, color = sample)) +
  geom_point(size = 3) + geom_text_repel(aes(label = id), seed = 1, max.overlaps = Inf, size = 3) +
  coord_equal() + labs(title = 'Haplotype PCoA',
    subtitle = str_glue('Distance = 1 - Jaccard; negative inertia {round(100 * negative_fraction, 3)}%'),
    x = axis_label(1), y = axis_label(2), color = 'Individual')
save_plot('haplotype_pcoa', ordination)

tree_status <- 'Not computed: neighbour joining needs at least three haplotypes.'
tree_plot <- ggplot() + annotate('text', x = 0, y = 0, label = tree_status) + theme_void()
if (n >= 3 && max(distance) > tol) {
  tree <- nj(as.dist(distance))
  write.tree(tree, file = file.path(out, 'haplotype_nj.nwk'))
  # Topological layout avoids distorting negative branches or implying an ancestral root.
  # Exact inferred lengths, including negatives, remain in Newick and the edge table.
  edges <- as_tibble(tree$edge, .name_repair = ~ c('parent', 'child')) |>
    mutate(length = tree$edge.length)
  tips <- seq_len(n)
  total <- n + tree$Nnode
  xx <- yy <- rep(NA_real_, total)
  root <- setdiff(edges$parent, edges$child) |> unique()
  if (length(root) != 1L) stop('Invalid NJ layout')
  tip_order <- 0
  position <- function(node, depth = 0) {
    xx[[node]] <<- depth
    children <- edges |> filter(parent == node) |> pull(child)
    if (!length(children)) {
      tip_order <<- tip_order + 1
      yy[[node]] <<- tip_order
    } else {
      walk(children, ~ position(.x, depth + 1))
      yy[[node]] <<- mean(yy[children])
    }
  }
  position(root)
  plotted <- edges |> mutate(x = xx[parent], y = yy[parent], xend = xx[child], yend = yy[child])
  write_tsv(edges, file.path(out, 'nj_edges.tsv'))
  labels <- tibble(id = tree$tip.label, x = xx[tips], y = yy[tips]) |>
    left_join(identity |> select(id, sample), by = 'id')
  tree_status <- str_glue('Unrooted NJ; display attachment is arbitrary; {sum(edges$length < -tol)} negative edges retained in Newick.')
  tree_plot <- ggplot(plotted) +
    geom_segment(aes(x = x, y = y, xend = x, yend = yend)) +
    geom_segment(aes(x = x, y = yend, xend = xend, yend = yend)) +
    geom_text(data = labels, aes(x = x + 0.1, y = y, label = id, color = sample), hjust = 0, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.65))) +
    theme_void() + labs(title = 'Descriptive neighbour-joining topology',
      subtitle = 'Branch lengths not drawn to scale; exact lengths supplied in Newick', color = 'Individual')
} else if (n >= 3) {
  tree_status <- 'Not computed: all distances are zero; topology is unresolved.'
  tree_plot <- ggplot() + annotate('text', x = 0, y = 0, label = tree_status) + theme_void()
}
save_plot('haplotype_comparison', (ordination / tree_plot) + plot_layout(guides = 'collect'), 11, 12)
write_lines(c('# Haplotype comparison', '',
  'CLIP distinct-sequence bp-weighted Jaccard similarity from Panacus. Distance = 1 - similarity. These views are descriptive, not population tests or species-tree inference. Clipping, assembly fragmentation and missing sequence can affect them.', '',
  '<details><summary>Similarity, ordination and tree</summary>', '',
  '![Similarity](haplotype_similarity.png)', '', '![PCoA and NJ](haplotype_comparison.png)', '',
  str_glue('PCoA retains {length(positive)} positive axes; negative inertia fraction = {negative_fraction}. Display percentages use positive inertia. No distance correction applied.'), '',
  tree_status, '',
  '[Distance matrix](haplotype_distance.tsv) · [PCoA coordinates](pcoa_coordinates.tsv) · [Eigenvalues](pcoa_eigenvalues.tsv)', '',
  if (file.exists(file.path(out, 'haplotype_nj.nwk'))) '[NJ Newick](haplotype_nj.nwk)' else '',
  '', '</details>'), file.path(out, 'comparison_report.md'))
packages <- c('tidyverse', 'ggplot2', 'patchwork', 'ape', 'ggrepel')
tibble(tool = c('R', packages), version = c(as.character(getRversion()), map_chr(packages, ~ as.character(packageVersion(.x))))) |>
  mutate(process = 'PANGENOME_COMPARISON_PLOTS', .before = 1) |>
  write_tsv(file.path(out, 'plot_versions.tsv'))
