/*
========================================================================================
    PANGENOME REPORT FIGURES MODULE  (workstream D)
========================================================================================
    Repo location: modules/pangenome_plots.nf

    Renders the pangenome report figures from the panacus coverage histogram and the
    variant catalog, via r_scripts/pangenome_plots.R. Growth/core curves, the Heaps'-law
    fit, and the confidence band are computed from the coverage histogram (rarefaction
    formulas); the SV size spectrum and variant-class bar come from the catalog.

    Reuses the pipeline R stack (ggplot2) via the 'pairwise_alignment' label.

    Input : tuple(taxid, hist_tsv, sv_sizes_tsv, variant_summary_tsv), plots_script
    Output: figures (PDFs) / growth_fit (machine-readable gamma / open-closed / sizes)
========================================================================================
*/

process PANGENOME_PLOTS {
    tag "${taxid}"
    label 'pairwise_alignment'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(hist), path(sv_sizes), path(variant_summary), path(hap_private)
    path(plots_script)

    output:
    tuple val(taxid), path("${taxid}.*.png"),          emit: figures, optional: true
    tuple val(taxid), path("${taxid}.growth_fit.tsv"), emit: growth_fit

    script:
    """
    Rscript ${plots_script} ${hist} ${sv_sizes} ${variant_summary} ${taxid} . \\
        hap_private=${hap_private} \\
        core=${params.pangenome_tier_core} \\
        softcore=${params.pangenome_tier_softcore} \\
        shell=${params.pangenome_tier_shell}
    """

    stub:
    """
    : > ${taxid}.growth_curves.png
    printf 'metric\\tvalue\\n' > ${taxid}.growth_fit.tsv
    """
}
