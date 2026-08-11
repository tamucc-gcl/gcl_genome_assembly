/*
========================================================================================
    PANGENOME PROGRESSIVE PLOT MODULE  (workstream H)
========================================================================================
    Repo location: modules/pangenome_progressive_plot.nf

    Renders the empirical progressive-growth curve from PANGENOME_PROGRESSIVE's table, via
    r_scripts/pangenome_progressive.R. Reuses 'pairwise_alignment' (ggplot2 only, no new
    package -> no cache bust of that env).

    Input : tuple(taxid, progressive_growth_tsv), progressive_script
    Output: png (progressive growth curve)
========================================================================================
*/

process PANGENOME_PROGRESSIVE_PLOT {
    tag "taxid_${taxid}"
    label 'pairwise_alignment'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(growth_tsv)
    path(progressive_script)

    output:
    tuple val(taxid), path("${taxid}.progressive_growth.png"), emit: png

    script:
    """
    Rscript ${progressive_script} ${growth_tsv} ${taxid} .
    """

    stub:
    """
    : > ${taxid}.progressive_growth.png
    """
}
