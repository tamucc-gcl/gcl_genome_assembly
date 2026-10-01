process PANGENOME_SHARING {
    tag "${taxid}"
    label 'panacus'
    conda 'bioconda::panacus=0.5.2 conda-forge::python=3.11'
    publishDir "${params.outdir}/pangenome/${taxid}/sharing", mode: params.publish_dir_mode, pattern: 'sharing/*', saveAs: { it.replaceFirst('^sharing/', '') }
    input:
    tuple val(taxid), path(gfa), path(identities)
    path(adapter)
    output:
    tuple val(taxid), path('sharing/sharing_report.md'), emit: report
    tuple val(taxid), path('sharing/sharing_summary.tsv'), emit: summary
    path('sharing/versions.tsv'), emit: versions
    tuple val(taxid), path('sharing/*'), emit: files
    script:
    """
    set -euo pipefail
    python3 ${adapter} prepare --gfa ${gfa} --ledger ${identities} --out sharing
    for unit in haplotype individual; do
        panacus -t ${task.cpus} hist -c bp --groupby sharing/\${unit}.groups.tsv ${gfa} > sharing/\${unit}.hist.tsv
    done
    python3 ${adapter} summarize --directory sharing \\
        --core ${params.pangenome_tier_core} --softcore ${params.pangenome_tier_softcore} --shell ${params.pangenome_tier_shell}
    for unit in haplotype individual; do
        panacus -t ${task.cpus} growth -c bp -l '${params.pangenome_growth_coverage}' -q '${params.pangenome_growth_quorum}' \\
            --groupby sharing/\${unit}.groups.tsv ${gfa} > sharing/\${unit}.growth.tsv
        test -s sharing/\${unit}.growth.tsv
    done
    python3 ${adapter} validate-growth --directory sharing
    printf 'process\\ttool\\tversion\\n' > sharing/versions.tsv
    printf '%s\\tpanacus\\t%s\\n' '${task.process}' "\$(panacus --version | sed -n '1p')" >> sharing/versions.tsv
    """
}
