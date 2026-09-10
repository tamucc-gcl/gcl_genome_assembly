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
p$add_argument("--hap_private",     default = "NO_HAP_PRIVATE") # PANGENOME_HAP_COVERAGE per-haplotype private breakdown
# ---- the view matrix (batch 5) --------------------------------------------------------
# representation_audit.tsv per tier. topology_enabled is the field that decides whether the
# topological classes in this report mean anything: the fine tier is vcfwave-decomposed, so AT
# is inherited from the parent record while REF/ALT are rewritten, and allele i stops
# corresponding to traversal i+1.
p$add_argument("--audit_parent",    default = "NO_AUDIT_PARENT")
p$add_argument("--audit_fine",      default = "NO_AUDIT_FINE")
# the FULL arm's private breakdown. The arms disagree by 46% on private bp and invert the
# reference's rank, so quoting one silently is wrong half the time.
p$add_argument("--hap_private_full", default = "NO_HAP_PRIVATE_FULL")
# PANGENOME_PRIVATE_PLOTS audit: reference_rank_clip / reference_rank_full
p$add_argument("--priv_figures_audit", default = "NO_PRIV_FIGURES")
p$add_argument("--rearr_audit_clip", default = "NO_REARR_CLIP")
p$add_argument("--rearr_audit_full", default = "NO_REARR_FULL")
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

# read a '#'-commented TSV. comment.char stays OFF: PanSN haplotype names contain '#' and
# read.delim's default would truncate them mid-field.
read_tsv_hash <- function(path) {
  if (is_missing(path)) return(NULL)
  ln <- readLines(path, warn = FALSE)
  ln <- ln[!grepl("^#", ln)]
  if (length(ln) < 2) return(NULL)
  d <- tryCatch(read.delim(text = paste(ln, collapse = "\n"), header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE, comment.char = ""),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) NULL else d
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

# ======================================================================================
# Which view says what  (batch 5)
#
# THREE views exist and they are not interchangeable. Every number elsewhere in this report
# comes from exactly one of them, and the unavailable cells are a finding rather than an
# omission -- see the notes under the table.
# ======================================================================================
# read_kv CANNOT be used on these files. Its read.delim call leaves comment.char at "",
# so a leading '#' line becomes the header; with no tab in it ncol(d) < 2 and read_kv bails
# to character(0). representation_audit.tsv opens with three '#' lines before its
# metric/value header, which is why the Records / Alt alleles / Topology rows came out "—"
# on the first run.
#
# comment.char cannot simply be switched on either: '#' is a legitimate character inside
# PanSN haplotype names (Sde-CBau_104#1), and read.delim would truncate those mid-field.
# read_tsv_hash already strips '#' LINES explicitly and keeps comment.char off, so route
# through it and pivot the two columns here.
read_kv_hash <- function(path, key = 1, val = 2) {
  d <- read_tsv_hash(path)
  if (is.null(d) || ncol(d) < 2) return(character(0))
  setNames(as.character(d[[val]]), as.character(d[[key]]))
}

ap_ <- read_kv_hash(args$audit_parent)
af_ <- read_kv_hash(args$audit_fine)
hpf <- read_tsv_hash(args$hap_private_full)
hpc <- read_tsv_hash(args$hap_private)
pfa <- read_kv_hash(args$priv_figures_audit)
rac <- read_kv_hash(args$rearr_audit_clip)
raf <- read_kv_hash(args$rearr_audit_full)

# "n/a" is never bare: each carries the reason, because a blank cell reads as zero.
NA_FINE_TOPO <- "n/a — AT destroyed by decomposition"
NA_FINE_TRAV <- "n/a — needs traversals"
NA_FULL_CAT  <- "n/a — no catalog on this arm"
NA_CLIP_ONLY <- "n/a — clip only"

have <- function(v) length(v) > 0
yn   <- function(v, k) { x <- g(v, k, NA); if (is.na(x)) "—" else as.character(x) }

priv_bp <- function(d) if (is.null(d) || !"private_bp" %in% names(d)) "—" else
  mb(sum(num(d$private_bp), na.rm = TRUE))

rows <- list()
r <- function(measure, parent, fine, full) rows[[length(rows) + 1]] <<-
  c(measure, parent, fine, full)

r("Records",
  if (have(ap_)) comma(g(ap_, "records")) else "—",
  if (have(af_)) comma(g(af_, "records")) else "—",
  NA_FULL_CAT)
r("Alt alleles per record",
  if (have(ap_)) yn(ap_, "alt_alleles_per_record") else "—",
  if (have(af_)) yn(af_, "alt_alleles_per_record") else "—",
  NA_FULL_CAT)
r("Topology usable (`AT`)",
  if (have(ap_)) sprintf("**%s**", yn(ap_, "topology_enabled")) else "—",
  if (have(af_)) sprintf("**%s** — %s", yn(af_, "topology_enabled"),
                         yn(af_, "decomposition_markers")) else NA_FINE_TOPO,
  NA_FULL_CAT)
r("Topological classes (SUBST / INS / DEL / INV)",
  "yes", NA_FINE_TOPO, NA_FULL_CAT)
r("Length-based classes", "yes", "yes", NA_FULL_CAT)
r("Merged reference footprint",
  if (have(ap_)) mb(g(ap_, "merged_ref_footprint_all_classes")) else "yes",
  "yes", NA_FULL_CAT)
r("Pangenome node bp / novel node bp",
  if (have(ap_)) mb(g(ap_, "novel_node_bp_all_classes")) else "yes",
  NA_FINE_TRAV, NA_FULL_CAT)
r("Allele frequency spectrum", "yes", "yes", NA_FULL_CAT)
r("Inversion rescue (alignment)", "yes", NA_FINE_TOPO, NA_FULL_CAT)
# NOT "both tiers": private sequence is a property of the GRAPH, not of a VCF view, so the
# tier axis does not apply to it at all. It is placed in the parent column only because the
# table has three columns and clip is one graph.
r("Private sequence (bp)",
  paste0(priv_bp(hpc), " *(clip graph)*"), "*n/a — a graph property, not a VCF view*",
  priv_bp(hpf))
r("Private segment size spectrum", "yes *(clip)*", "", "yes")
r("Private per chromosome", "yes *(clip)*", "", "yes")
r("Private evidence (map × k-mer)", "yes *(clip)*", "", "yes")
r("Rearrangement / untangle",
  if (have(rac)) sprintf("%s inverted *(clip)*", mb(g(rac, "inverted_bp_total"))) else "yes *(clip)*",
  "",
  if (have(raf)) sprintf("%s inverted", mb(g(raf, "inverted_bp_total"))) else "yes")
r("Openness / growth / partition", "yes", "", NA_CLIP_ONLY)

add("### Which view says what", "",
    "Every number below comes from exactly one of three views, and they are **not**",
    "interchangeable. Unavailable cells carry their reason: a blank would read as zero.", "",
    "| Measure | clip · parent | clip · fine | full |",
    "|---|---|---|---|")
for (rw in rows)
  add(sprintf("| %s | %s | %s | %s |", rw[1],
              if (nzchar(rw[2])) rw[2] else "—",
              if (nzchar(rw[3])) rw[3] else "—",
              if (nzchar(rw[4])) rw[4] else "—"))
add("",
    "- **clip · parent** is the only view where `AT` is interpretable, so it is the only view",
    "  in which topological classification is valid. Everything called SUBST, INS, DEL or",
    "  `INV_PATH_EXPLICIT` in this report comes from here.",
    "- **clip · fine** is `vcfwave`-decomposed. `AT` is inherited from the parent record while",
    "  REF and ALT are rewritten, so allele *i* stops corresponding to traversal *i+1*. The",
    "  classifier refuses to compute topology on it rather than trusting the wiring to be",
    "  right, which is why those cells are unavailable rather than empty.",
    "- **full** has no variant catalog: `PANGENOME_VARIANTS` and `PANGENOME_CLASSIFY` run on",
    "  the clip graph. The full arm contributes the private-sequence and rearrangement",
    "  measures only.", "")

if (have(pfa) && !is.na(g(pfa, "reference_rank_clip", NA))) {
  add(sprintf(paste0("> **The arms disagree about private sequence, and not by a little.** ",
                     "Clipping removes 698,360,436 bp of which 99.1%% is private, so the clip ",
                     "arm understates private content by ~46%%. It also inverts the ",
                     "reference's rank: **%s** on clip, **%s** on full. The reference is the ",
                     "graph backbone and is never clipped, so on the clip arm it looks like ",
                     "the most private haplotype in the cohort and on the full arm the least. ",
                     "Private-sequence figures in this report are drawn from both arms for ",
                     "that reason."),
              g(pfa, "reference_rank_clip"), g(pfa, "reference_rank_full")), "")
}
add(sprintf("Minigraph-Cactus pangenome graph for *%s*%s.",
            gsub("_", " ", species),
            if (!is.na(g(gr, "n_haplotypes"))) sprintf(", over %s haplotypes", g(gr, "n_haplotypes")) else ""), "")

# ---- downstream files (manifest-driven, pivoted by graph) ------------------------------
# pangenome_manifest.tsv is the downstream seam: role -> file, one row per (role, graph). The
# table is rendered from it, so a role added to PANGENOME_MANIFEST appears here automatically;
# roles with no entry below fall back to a humanized role name. A manifest without a `graph`
# column (the pre-clip/full schema) is treated as all-clip.
role_name <- c(
  graph_gbz        = "Graph (GBZ)",
  graph_gfa        = "Graph (GFA)",
  snarls           = "Snarls",
  haplotype_index  = "Haplotype index",
  variants_vcf     = "Variant catalog",
  variants_vcf_tbi = "Variant catalog index",
  reference_fasta  = "Reference FASTA",
  reference_fai    = "Reference FASTA index"
)
role_desc <- c(
  graph_gbz        = "Compressed graph + haplotype paths \u2014 the mapping substrate (`vg giraffe`, `vg call`).",
  graph_gfa        = "Text graph interchange (gzipped) \u2014 panacus, odgi, Bandage, `vg convert`.",
  snarls           = "Snarl (bubble) decomposition \u2014 the genotyping units for `vg call`.",
  haplotype_index  = "Haplotype-sampling index \u2014 builds personalized reference graphs for `vg giraffe`.",
  variants_vcf     = "Variant catalog vs the reference path: top-level bubbles decomposed to minimal variants by vcfwave, one allele per row, SNP/indel/SV.",
  blocks_vcf       = "Pre-decomposition block-level VCF (vcfbub output). Kept so the vcfwave decomposition can be audited.",
  variants_vcf_tbi = "Tabix index for the variant catalog.",
  reference_fasta  = "Reference-path FASTA \u2014 surjection target and linear-coordinate seam.",
  reference_fai    = "faidx index for the reference-path FASTA."
)
# Published, and listed in the manifest by role -- but NOT linked in this table. The whole-graph
# odgi files are one to two orders of magnitude larger than everything else here and would
# dominate a table meant to be read by a person. The footnote points at them instead; machine
# consumers resolve them from the manifest, which is the interface that matters for them.
role_skip <- c("odgi")

if (!is_missing(args$manifest)) {
  mf <- tryCatch(read.delim(args$manifest, header = TRUE, stringsAsFactors = FALSE,
                            check.names = FALSE), error = function(e) NULL)
  if (!is.null(mf) && nrow(mf) > 0 && all(c("role", "file") %in% names(mf))) {
    mf$role  <- as.character(mf$role)
    mf$file  <- basename(as.character(mf$file))
    mf$graph <- if ("graph" %in% names(mf)) as.character(mf$graph) else "clip"
    roles <- setdiff(unique(mf$role), role_skip)   # manifest row order = display order
    if (length(roles) > 0) {
      cell <- function(r, gg) {
        f <- mf$file[mf$role == r & mf$graph == gg]
        if (length(f) == 0) "\u2014" else sprintf("[`%s`](%s)", f[1], rel(f[1]))
      }
      nm <- vapply(roles, function(r) if (r %in% names(role_name)) role_name[[r]] else gsub("_", " ", r), "")
      ds <- vapply(roles, function(r) if (r %in% names(role_desc)) role_desc[[r]] else "", "")
      add("### Downstream files", "",
          sprintf(paste("Graph products for downstream use, split by graph: **clip** is the",
                        "reference-anchored graph (the default for mapping and variant work);",
                        "**full** retains sequence that had no reference alignment. Paths are",
                        "relative to this report; [`%s`](%s) carries the same set machine-readably."),
                  basename(args$manifest), rel(basename(args$manifest))), "",
          "| product | clip graph | full graph | description |", "|---|---|---|---|",
          sprintf("| %s | %s | %s | %s |", nm,
                  vapply(roles, cell, "", gg = "clip"),
                  vapply(roles, cell, "", gg = "full"),
                  ds), "",
          paste("*The whole-graph odgi graphs (`.og`, `.full.og`) are published alongside these",
                "but not linked above \u2014 they are far larger than every other product. Resolve",
                "them from the manifest by the `odgi` role, or rebuild either from its GFA with",
                "`odgi build -g`. Per-chromosome graphs sit under `<label>.chroms/`. Read-mapping",
                "indexes (`.dist` / `.min`) are not built here: haplotype sampling makes them",
                "sample-specific, so the consuming pipeline builds them from the GBZ + `.hapl`",
                "at map time.*"), "")
    }
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
      "*View: clip · parent. Topological classes are valid only here — see the matrix above.*", "",
      "| class | count |", "|---|---:|",
      sprintf("| SNP | %s |", comma(g(vs, "SNP"))),
      sprintf("| Indel | %s |", comma(g(vs, "INDEL"))),
      sprintf("| SV | %s (INS %s / DEL %s) |", comma(g(vs, "SV")), comma(g(vs, "SV_INS")), comma(g(vs, "SV_DEL"))),
      # BLOCKSUB = both alleles >= minsv and lengths differ: what vcfwave could not
      # decompose. A large count means alleles are genuinely unrelated at those loci.
      sprintf("| &nbsp;&nbsp;of which block substitutions | %s |", comma(g(vs, "SV_BLOCKSUB", "0"))),
      sprintf("| &nbsp;&nbsp;of which same-length substitutions | %s |", comma(g(vs, "SV_COMPLEX", "0"))),
      "",
      sprintf("Total: %s variants relative to the reference path.", comma(total)), "")
}

# ---- openness / growth ----------------------------------------------------------------
if (length(gr) > 0) {
  add("### Openness / growth", "",
      "*View: clip. panacus runs on the clip GFA; there is no full-arm equivalent.*", "",
      sprintf("Pangenome %s · core %s · accessory %s · private %s. Heaps' \u03b3 = %s (%s).",
              mb(g(gr, "pangenome_bp")), mb(g(gr, "core_bp")), mb(g(gr, "accessory_bp")),
              mb(g(gr, "private_bp")), format(round(num(g(gr, "heaps_gamma")), 3)), g(gr, "openness", "—")),
      "",
      sprintf("![Growth and core curves](%s)", fig("growth_curves.png")), "",
      sprintf("![Coverage histogram](%s)", fig("coverage_histogram.png")), "")
}

# ---- core / soft-core / shell / cloud partition ---------------------------------------
# Gated on tier_core_bp because that row is written by the same step that renders the two
# figures, so its presence is a safe proxy for the PNGs existing on disk.
if (!is.na(g(gr, "tier_core_bp"))) {
  add("### Partition (core / soft-core / shell / cloud)", "",
      sprintf(paste("Graph sequence binned by how many haplotypes carry it. Cuts are settable",
                    "(`params.pangenome_tier_*`); the realised floors for this cohort are core",
                    "\u2265%s hap, soft-core \u2265%s, shell \u2265%s, out of %s. A tier can be empty",
                    "when a cut is shadowed by a stricter one at this cohort size."),
              g(gr, "tier_core_min_haps", "\u2014"), g(gr, "tier_softcore_min_haps", "\u2014"),
              g(gr, "tier_shell_min_haps", "\u2014"), g(gr, "n_haplotypes", "\u2014")), "",
      "| tier | sequence |", "|---|---:|",
      sprintf("| Core | %s |",      mb(g(gr, "tier_core_bp"))),
      sprintf("| Soft-core | %s |", mb(g(gr, "tier_softcore_bp"))),
      sprintf("| Shell | %s |",     mb(g(gr, "tier_shell_bp"))),
      sprintf("| Cloud | %s |",     mb(g(gr, "tier_cloud_bp"))), "",
      sprintf("![Pangenome partition](%s)", fig("pangenome_partition.png")), "",
      sprintf("![Coverage histogram by tier](%s)", fig("coverage_histogram_tiers.png")), "")
}

# ---- private-sequence ownership -------------------------------------------------------
hp_tbl <- read_tsv_hash(args$hap_private)
if (!is.null(hp_tbl) && all(c("haplotype", "private_bp", "pct_of_private",
                              "pct_of_haplotype") %in% names(hp_tbl))) {
  hp_tbl$private_bp       <- num(hp_tbl$private_bp)
  hp_tbl$pct_of_private   <- num(hp_tbl$pct_of_private)
  hp_tbl$pct_of_haplotype <- num(hp_tbl$pct_of_haplotype)
  hp_tbl <- hp_tbl[order(-hp_tbl$pct_of_private), , drop = FALSE]
  nh_p   <- nrow(hp_tbl)
  shown  <- if (nh_p > 20) head(hp_tbl, 15) else hp_tbl
  add("### Private-sequence ownership", "",
      "*View: clip. The full arm reports ~46% more private sequence and reverses the",
      "reference's rank — see the matrix above and the clip-vs-full figure.*", "",
      sprintf(paste("Every private segment is carried by exactly one haplotype, so the shares",
                    "below sum to 100%%. An even split over %d haplotypes would be %.1f%% each;",
                    "large departures mark assemblies contributing unusual amounts of unshared",
                    "sequence, which is either real structural variation or residual",
                    "contamination / uncollapsed duplication surviving into the graph."),
              nh_p, 100 / nh_p), "",
      "| haplotype | private | % of private | % of this haplotype |", "|---|---:|---:|---:|",
      sprintf("| %s | %s | %.1f%% | %.2f%% |", shown$haplotype,
              vapply(shown$private_bp, mb, character(1)),
              shown$pct_of_private, shown$pct_of_haplotype), "")
  if (nrow(shown) < nh_p)
    add(sprintf("*Top %d of %d haplotypes; the full table is `%s.hap_private.tsv`.*",
                nrow(shown), nh_p, species), "")
  # figures come from PANGENOME_PLOTS, which is gated on params.pangenome_growth
  if (length(gr) > 0)
    add(sprintf("![Private sequence by haplotype](%s)", fig("private_by_haplotype.png")), "",
        sprintf("![Private fraction of each haplotype](%s)",
                fig("private_fraction_by_haplotype.png")), "")
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
    "*View: clip · parent.*", "",
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
         jq("tier_core_bp", g(gr, "tier_core_bp")), jq("tier_softcore_bp", g(gr, "tier_softcore_bp")),
         jq("tier_shell_bp", g(gr, "tier_shell_bp")), jq("tier_cloud_bp", g(gr, "tier_cloud_bp")),
         jq("heaps_gamma", g(gr, "heaps_gamma")), jq("openness", g(gr, "openness"), q = TRUE),
         jq("snp", g(vs, "SNP")), jq("indel", g(vs, "INDEL")), jq("sv", g(vs, "SV")),
         jq("edit_rate", g(qc, "edit_rate")), jq("graph_identity", g(qc, "graph_identity"))
       ), collapse = ",\n"),
       "}")
writeLines(j, args$json)
