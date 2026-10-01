process PANGENOME_SHARING_PLOTS {
    tag "${taxid}"
    cpus 1
    memory '2 GB'
    time '30m'
    conda 'conda-forge::r-base=4.3.3'
    publishDir "${params.outdir}/pangenome/${taxid}/figures", mode: params.publish_dir_mode,
        pattern: 'figures/*', saveAs: { it.replaceFirst('^figures/', '') }
    input:
    tuple val(taxid), path(tables, stageAs: 'sharing/*')
    path(plotter)
    output:
    tuple val(taxid), path('figures/sharing_figures.md'), emit: report
    tuple val(taxid), path('figures/*.png'), emit: figures
    path('figures/*.pdf'), emit: pdfs
    path('figures/versions.tsv'), emit: versions
    script:
    """
    set -euo pipefail
    Rscript ${plotter} sharing figures
    """
}
