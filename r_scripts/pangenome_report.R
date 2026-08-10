#!/usr/bin/env Rscript
# ======================================================================================
# pangenome_report.R  (workstream F)
# Repo location: r_scripts/pangenome_report.R
#
# Assembles a SELF-CONTAINED pangenome report section (markdown) + a machine-readable
# stats JSON from the pangenome stats tables. The markdown is written in the same style as
# generate_summary_report.R (a character vector) so it can either stand alone or be appended
# as a knitr-free child section by the main report (via readLines).
#
# Image paths are written relative to the pipeline outdir root (pangenome/<species>/<png>)
# so they resolve when the fragment is included in the top-level assembly_report.md.
#
# Usage:
#   Rscript pangenome_report.R --qc_metrics X --growth_fit X --variant_summary X \
#       --graph_stats X --species NAME --output pangenome_report.md --json pangenome_stats.json
# Any input may be a sentinel (NO_*) or missing/empty; that part is simply skipped.
# ======================================================================================

suppressPackageStartupMessages({ library(argparse) })

p <- ArgumentParser()
p$add_argument("--qc_metrics",      default = "NO_QC")
p$add_argument("--growth_fit",      default = "NO_GROWTH")
p$add_argument("--variant_summary", default = "NO_VARIANTS")
p$add_argument("--graph_stats",     default = "NO_GRAPH")   # odgi stats -S text (length/nodes/edges/paths/steps)
p$add_argument("--popstruct",       default = "NO_POPSTRUCT") # sentinel/indicator; PCA+NJ PNGs sit next to the other figures
p$add_argument("--species",         default = "pangenome")
p$add_argument("--output",          default = "pangenome_report.md")
p$add_argument("--json",            default = "pangenome_stats.json")
args <- p$parse_args()

is_missing <- function(x) is.null(x) || is.na(x) || grepl("^NO_", basename(x)) || !file.exists(x) || file.size(x) == 0

# read a two-column key/value TSV -> named character vector
read_kv <- function(path, key = 1, val = 2) {
  if (is_missing(path)) return(character(0))
  d <- tryCatch(read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || ncol(d) < 2) return(character(0))
  setNames(as.character(d[[val]]), as.character(d[[key]]))
}
g   <- function(v, k, d = NA) if (k %in% names(v)) v[[k]] else d
num <- function(x) suppressWarnings(as.numeric(x))
mb  <- function(x) { x <- num(x); if (is.na(x)) "—" else paste0(format(round(x/1e6, 1), big.mark = ","), " Mb") }
gb  <- function(x) { x <- num(x); if (is.na(x)) "—" else paste0(format(round(x/1e9, 3), nsmall = 3), " Gb") }
comma <- function(x) { x <- num(x); if (is.na(x)) "—" else format(round(x), big.mark = ",", scientific = FALSE) }

qc  <- read_kv(args$qc_metrics)
gr  <- read_kv(args$growth_fit)
vs  <- read_kv(args$variant_summary)   # class -> count

# graph stats: parse the odgi stats -S text (header '#length nodes edges paths steps' + values row)
graph <- list(length = NA, nodes = NA, edges = NA, paths = NA, steps = NA)
if (!is_missing(args$graph_stats)) {
  ln <- readLines(args$graph_stats, warn = FALSE)
  hdr <- grep("^#length", ln)
  if (length(hdr) >= 1 && length(ln) > hdr[1]) {
    vals <- strsplit(trimws(ln[hdr[1] + 1]), "\\s+")[[1]]
    if (length(vals) >= 5) graph <- list(length = vals[1], nodes = vals[2], edges = vals[3],
                                          paths = vals[4], steps = vals[5])
  }
}

species <- args$species
fig <- function(suffix) sprintf("pangenome/%s/%s.%s", species, species, suffix)
md <- character(0)
add <- function(...) md <<- c(md, ...)

nhap <- g(gr, "n_haplotypes", g(qc, "n_alignments", NA))
add("## Pangenome", "")
add(sprintf("Minigraph-Cactus pangenome graph for *%s*%s.",
            gsub("_", " ", species),
            if (!is.na(g(gr, "n_haplotypes"))) sprintf(", over %s haplotypes", g(gr, "n_haplotypes")) else ""), "")

# ---- graph ----------------------------------------------------------------------------
add("### Graph", "",
    "| property | value |", "|---|---:|",
    sprintf("| Haplotypes | %s |", g(gr, "n_haplotypes", "—")),
    sprintf("| Total length | %s |", gb(graph$length)),
    sprintf("| Nodes | %s |", comma(graph$nodes)),
    sprintf("| Edges | %s |", comma(graph$edges)),
    sprintf("| Paths | %s |", comma(graph$paths)),
    sprintf("| Acyclic | %s |", g(qc, "is_acyclic", "—")), "")

# ---- variant catalog ------------------------------------------------------------------
if (length(vs) > 0) {
  total <- sum(num(vs[c("SNP","INDEL","SV")]), na.rm = TRUE)
  add("### Variant catalog", "",
      "| class | count |", "|---|---:|",
      sprintf("| SNP | %s |", comma(g(vs, "SNP"))),
      sprintf("| Indel | %s |", comma(g(vs, "INDEL"))),
      sprintf("| SV | %s (INS %s / DEL %s) |", comma(g(vs, "SV")), comma(g(vs, "SV_INS")), comma(g(vs, "SV_DEL"))),
      "",
      sprintf("Total: %s variants relative to the reference path.", comma(total)), "")
}

# ---- openness / growth ----------------------------------------------------------------
if (length(gr) > 0) {
  add("### Openness / growth", "",
      sprintf("Pangenome %s · core %s · accessory %s · private %s. Heaps' \u03b3 = %s (%s).",
              mb(g(gr, "pangenome_bp")), mb(g(gr, "core_bp")), mb(g(gr, "accessory_bp")),
              mb(g(gr, "private_bp")), format(round(num(g(gr, "heaps_gamma")), 3)), g(gr, "openness", "—")),
      "",
      sprintf("![Growth and core curves](%s)", fig("growth_curves.png")), "",
      sprintf("![Coverage histogram](%s)", fig("coverage_histogram.png")), "")
}

# ---- structural variants figure -------------------------------------------------------
add("### Structural variants", "",
    sprintf("![SV size spectrum](%s)", fig("sv_size_histogram.png")), "")

# ---- population structure (PCoA + NJ tree, per haplotype and per individual) ----------
if (!is_missing(args$popstruct)) {
  add("### Population structure", "",
      sprintf(paste("Ordination (PCoA) and neighbour-joining trees from graph-similarity",
                    "distances (shared node content \u2014 SNPs, indels and SVs), at two levels:",
                    "per haploid assembly (%s haplotypes, incl. the reference) and per diploid",
                    "individual (haplotypes aggregated)."),
              g(gr, "n_haplotypes", "sampled")),
      "",
      "**Per haplotype**", "",
      sprintf("![Haplotype PCoA](%s)", fig("pca_haplotype.png")), "",
      sprintf("![Haplotype NJ tree](%s)", fig("njtree_haplotype.png")), "",
      "**Per individual**", "",
      sprintf("![Individual PCoA](%s)", fig("pca_individual.png")), "",
      sprintf("![Individual NJ tree](%s)", fig("njtree_individual.png")), "")
}

# ---- graph quality --------------------------------------------------------------------
if (length(qc) > 0) {
  add("### Graph quality", "",
      "| metric | value |", "|---|---:|",
      sprintf("| Re-alignment identity | %s (edit rate %s) |", g(qc, "graph_identity", "—"), g(qc, "edit_rate", "—")),
      sprintf("| Realigned | %s over %s alignments |", gb(g(qc, "realigned_bp")), comma(g(qc, "n_alignments"))),
      sprintf("| Mean node degree | %s (max %s) |", g(qc, "avg_node_degree", "—"), g(qc, "max_node_degree", "—")),
      sprintf("| Mean links length | %s bp |", g(qc, "mean_links_length", "—")), "",
      paste("*Re-alignment identity is minigraph-GAF level (coarse); per-haplotype assembly",
            "BUSCO/QV remain the authoritative per-assembly completeness/accuracy metrics.*"), "")
}

writeLines(md, args$output)

# ---- machine-readable stats JSON ------------------------------------------------------
jq  <- function(k, v, q = FALSE) {
  if (is.na(v) || v == "—") return(sprintf('  "%s": null', k))
  if (q) sprintf('  "%s": "%s"', k, v) else sprintf('  "%s": %s', k, v)
}
j <- c("{",
       paste(c(
         jq("species", species, q = TRUE),
         jq("n_haplotypes", g(gr, "n_haplotypes")),
         jq("graph_length_bp", graph$length),
         jq("nodes", graph$nodes), jq("edges", graph$edges), jq("paths", graph$paths),
         jq("pangenome_bp", g(gr, "pangenome_bp")), jq("core_bp", g(gr, "core_bp")),
         jq("accessory_bp", g(gr, "accessory_bp")), jq("private_bp", g(gr, "private_bp")),
         jq("heaps_gamma", g(gr, "heaps_gamma")), jq("openness", g(gr, "openness"), q = TRUE),
         jq("snp", g(vs, "SNP")), jq("indel", g(vs, "INDEL")), jq("sv", g(vs, "SV")),
         jq("edit_rate", g(qc, "edit_rate")), jq("graph_identity", g(qc, "graph_identity"))
       ), collapse = ",\n"),
       "}")
writeLines(j, args$json)
