process PANGENOME_SIMILARITY {
    tag "${taxid}"
    label 'panacus'
    publishDir "${params.outdir}/pangenome/${taxid}/comparison", mode: params.publish_dir_mode
    input:
    tuple val(taxid), path(gfa), path(identities), path(groups)
    output:
    tuple val(taxid), path('haplotype_similarity.tsv'), path(identities), emit: matrix
    path('versions.tsv'), emit: versions
    script:
    """
    set -euo pipefail
    panacus -t ${task.cpus} similarity -c bp --groupby ${groups} ${gfa} > haplotype_similarity.tsv
    test -s haplotype_similarity.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    printf '%s\\tpanacus\\t%s\\n' '${task.process}' "\$(panacus --version | sed -n '1p')" >> versions.tsv
    """
}

process PANGENOME_COMPARISON_PLOTS {
    tag "${taxid}"
    cpus 1
    memory '4 GB'
    time '1h'
    conda 'conda-forge::r-base=4.3.3 conda-forge::r-tidyverse=2.0.0 conda-forge::r-ggplot2=3.5.1 conda-forge::r-patchwork=1.3.0 conda-forge::r-ape=5.8 conda-forge::r-ggrepel=0.9.6'
    publishDir "${params.outdir}/pangenome/${taxid}/comparison", mode: params.publish_dir_mode,
        pattern: 'comparison/*', saveAs: { it.replaceFirst('^comparison/', '') }
    input:
    tuple val(taxid), path(similarity), path(identities)
    path(plotter)
    output:
    tuple val(taxid), path('comparison/comparison_report.md'), emit: report
    tuple val(taxid), path('comparison/*'), emit: files
    path('comparison/plot_versions.tsv'), emit: versions
    script:
    """
    set -euo pipefail
    Rscript ${plotter} ${similarity} ${identities} comparison
    """
}
