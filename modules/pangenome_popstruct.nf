/*
========================================================================================
    PANGENOME POPULATION-STRUCTURE PLOTS MODULE  (workstream D)
========================================================================================
    Repo location: modules/pangenome_popstruct.nf

    Renders the per-haplotype PCoA scatter + neighbour-joining tree from the odgi similarity
    matrix, via r_scripts/pangenome_popstruct.R. One point / one leaf per haplotype (incl.
    the reference). Always emits both PNGs (a labelled placeholder only if < 3 haplotypes).

    Dedicated 'pangenome_popstruct' env (r-base + r-ggplot2 + r-ape + r-ggrepel; cmdscale is
    base R) so adding these does not bust the shared 'pairwise_alignment' cache. Labels are
    shortened -- a species tag such as 'Sde-' is dropped when every unit shares it, and the
    reference individual's flat '<ind>_hapN#0' names collapse to '<ind>#N' so the reference is
    not split into two pseudo-individuals. PCoA labels are repelled (ggrepel) so they neither
    overlap nor run off the page; the NJ trees are midpoint-rooted rightwards phylograms with
    aligned tip labels (midpoint rooting is done with ape primitives -- no extra dependency).

    Input : tuple(taxid, similarity_tsv), popstruct_script
    Output: pca_png (PCoA) / nj_png
========================================================================================
*/

process PANGENOME_POPSTRUCT {
    tag "taxid_${taxid}"
    label 'pangenome_popstruct'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(similarity)
    path(popstruct_script)

    output:
    tuple val(taxid), path("${taxid}.pca_haplotype.png"), emit: pca_png    // report indicator
    tuple val(taxid), path("${taxid}.*.png"),             emit: figures    // all four (hap + individual)

    script:
    """
    Rscript ${popstruct_script} ${similarity} ${taxid} .
    """

    stub:
    """
    : > ${taxid}.pca_haplotype.png
    : > ${taxid}.njtree_haplotype.png
    : > ${taxid}.pca_individual.png
    : > ${taxid}.njtree_individual.png
    """
}
