process PANGENOME_NODE_COVERAGE {
    tag "${taxid}"
    label 'panacus'
    input:
    tuple val(taxid), path(gfa), path(groups)
    output:
    tuple val(taxid), path('node_coverage.tsv.gz'), emit: coverage
    path('versions.tsv'), emit: versions
    script:
    """
    set -euo pipefail
    panacus -t ${task.cpus} table -c bp --total --groupby ${groups} ${gfa} | gzip -1 > node_coverage.tsv.gz
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    printf '%s\\tpanacus\\t%s\\n' '${task.process}' "\$(panacus --version | sed -n '1p')" >> versions.tsv
    """
}

process PANGENOME_REGIONAL {
    tag "${taxid}"
    cpus 1
    memory '24 GB'
    time '24h'
    conda 'conda-forge::python=3.11 conda-forge::numpy=1.26.4'
    publishDir "${params.outdir}/pangenome/${taxid}/regional", mode: params.publish_dir_mode,
        pattern: 'regional/*', saveAs: { it.replaceFirst('^regional/', '') }
    input:
    tuple val(taxid), path(gfa), path(groups), path(histogram), path(coverage)
    path(adapter)
    path(regional_tests)
    output:
    tuple val(taxid), path('regional/regional_coverage.tsv'), emit: coverage
    tuple val(taxid), path('regional/regional_audit.json'), emit: audit
    tuple val(taxid), path('regional/regional_report.md'), emit: report
    tuple val(taxid), path('regional/regional_sharing.tsv'), emit: sharing
    script:
    """
    set -euo pipefail
    python3 ${regional_tests}
    python3 ${adapter} --gfa ${gfa} --totals ${coverage} --groups ${groups} --histogram ${histogram} --out regional \\
        --core ${params.pangenome_tier_core} --softcore ${params.pangenome_tier_softcore} --shell ${params.pangenome_tier_shell}
    """
}
