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
p$add_argument("--progressive",     default = "NO_PROGRESSIVE") # sentinel/indicator; progressive growth PNG (opt-in)
p$add_argument("--manifest",        default = "NO_MANIFEST")   # pangenome_manifest.tsv (role/file/label)
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
# Paths are written relative to the outdir ROOT, which is where the main report lives, so
# both figures and file links resolve once this fragment is embedded in assembly_report.md.
rel <- function(f)      sprintf("pangenome/%s/%s", species, f)
fig <- function(suffix) rel(sprintf("%s.%s", species, suffix))
md <- character(0)
add <- function(...) md <<- c(md, ...)

nhap <- g(gr, "n_haplotypes", g(qc, "n_alignments", NA))
add("## Pangenome", "")
add(sprintf("Minigraph-Cactus pangenome graph for *%s*%s.",
            gsub("_", " ", species),
            if (!is.na(g(gr, "n_haplotypes"))) sprintf(", over %s haplotypes", g(gr, "n_haplotypes")) else ""), "")

# ---- downstream files (manifest-driven; sits at the top of the section) ----------------
# pangenome_manifest.tsv is the frozen downstream seam (role -> file). Rendering the table
# from it means a new role added to PANGENOME_MANIFEST appears here automatically; roles
# without a description below still render (the role name is humanized as a fallback).
role_desc <- c(
  graph_gbz        = "Compressed graph + haplotype paths \u2014 the mapping substrate (`vg giraffe`, `vg call`).",
  graph_gfa        = "Text graph interchange (gzipped) \u2014 panacus, odgi, Bandage, `vg convert`.",
  odgi             = "odgi graph \u2014 path/node queries, similarity, sorting, 1D/2D layout.",
  snarls           = "Snarl (bubble) decomposition \u2014 the genotyping units for `vg call`.",
  haplotype_index  = "Haplotype-sampling index \u2014 builds personalized reference graphs for `vg giraffe`.",
  variants_vcf     = "Variant catalog vs the reference path: top-level bubbles, one allele per row, SNP/indel/SV.",
  variants_vcf_tbi = "Tabix index for the variant catalog.",
  reference_fasta  = "Reference-path FASTA \u2014 surjection target and linear-coordinate seam.",
  reference_fai    = "faidx index for the reference-path FASTA.",
  manifest         = "This table, machine-readable (`role` / `file` / `label`) \u2014 resolve files by role, not by name."
)

if (!is_missing(args$manifest)) {
  mf <- tryCatch(read.delim(args$manifest, header = TRUE, stringsAsFactors = FALSE,
                            check.names = FALSE), error = function(e) NULL)
  if (!is.null(mf) && nrow(mf) > 0 && all(c("role", "file") %in% names(mf))) {
    rows <- data.frame(role = as.character(mf$role),
                       file = basename(as.character(mf$file)),
                       stringsAsFactors = FALSE)
    if (!"manifest" %in% rows$role) {
      rows <- rbind(rows, data.frame(role = "manifest",
                                     file = basename(args$manifest),
                                     stringsAsFactors = FALSE))
    }
    desc <- ifelse(rows$role %in% names(role_desc),
                   role_desc[rows$role],
                   gsub("_", " ", rows$role))
    add("### Downstream files", "",
        paste("The graph products an external pipeline consumes. Paths are relative to this",
              "report; `pangenome_manifest.tsv` carries the same set machine-readably."),
        "",
        "| file | description |", "|---|---|",
        sprintf("| [`%s`](%s) | %s |", rows$file, rel(rows$file), desc), "",
        paste("*Read-mapping indexes (`.dist` / `.min`) are deliberately not built here \u2014",
              "haplotype sampling makes them sample-specific, so the consuming pipeline builds",
              "them from the GBZ + `.hapl` at map time.*"), "")
  }
}

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

# ---- progressive (incremental-construction) growth, if run (opt-in) --------------------
if (!is_missing(args$progressive)) {
  add("### Progressive growth (empirical)", "",
      paste("Graph size as each assembly is added (minigraph, reference-first) \u2014 the",
            "empirical counterpart to the analytic growth above."),
      "",
      sprintf("![Progressive growth](%s)", fig("progressive_growth.png")), "")
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
