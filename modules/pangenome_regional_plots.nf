process PANGENOME_REGIONAL_PLOTS {
    tag "${taxid}"
    cpus 1
    memory '4 GB'
    time '1h'
    conda 'conda-forge::r-base=4.3.3 conda-forge::r-tidyverse=2.0.0 conda-forge::r-ggplot2=3.5.1 conda-forge::r-patchwork=1.3.0'
    publishDir "${params.outdir}/pangenome/${taxid}/regional/figures", mode: params.publish_dir_mode,
        pattern: 'figures/*', saveAs: { it.replaceFirst('^figures/', '') }
    input:
    tuple val(taxid), path(sharing)
    path(plotter)
    output:
    tuple val(taxid), path('figures/regional_figures.md'), emit: report
    tuple val(taxid), path('figures/*'), emit: files
    path('figures/plot_versions.tsv'), emit: versions
    script:
    """
    set -euo pipefail
    Rscript ${plotter} ${sharing} figures
    """
}
