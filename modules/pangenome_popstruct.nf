/*
========================================================================================
    PANGENOME POPULATION-STRUCTURE PLOTS MODULE  (workstream D)
========================================================================================
    Repo location: modules/pangenome_popstruct.nf

    Renders the per-haplotype PCoA scatter + neighbour-joining tree from the odgi similarity
    matrix, via r_scripts/pangenome_popstruct.R. One point / one leaf per haplotype (incl.
    the reference). Always emits both PNGs (a labelled placeholder only if < 3 haplotypes).

    Reuses the pipeline R stack via 'pairwise_alignment' -- that env needs `ape` added
    (conda-forge::r-ape) alongside its existing ggplot2. cmdscale is base R.

    Input : tuple(taxid, similarity_tsv), popstruct_script
    Output: pca_png (PCoA) / nj_png
========================================================================================
*/

process PANGENOME_POPSTRUCT {
    tag "taxid_${taxid}"
    label 'pairwise_alignment'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(similarity)
    path(popstruct_script)

    output:
    tuple val(taxid), path("${taxid}.pca.png"),    emit: pca_png
    tuple val(taxid), path("${taxid}.njtree.png"), emit: nj_png

    script:
    """
    Rscript ${popstruct_script} ${similarity} ${taxid} .
    """

    stub:
    """
    : > ${taxid}.pca.png
    : > ${taxid}.njtree.png
    """
}
