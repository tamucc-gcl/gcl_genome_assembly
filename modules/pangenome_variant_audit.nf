process PANGENOME_VARIANT_AUDIT {
    tag "${taxid}"
    cpus 2
    memory '8 GB'
    time '8h'
    conda 'bioconda::bcftools=1.21 conda-forge::python=3.11'
    publishDir "${params.outdir}/pangenome/${taxid}/variants", mode: params.publish_dir_mode,
        pattern: 'variants/*', saveAs: { it.replaceFirst('^variants/', '') }
    input:
    tuple val(taxid), path(vcf), path(tbi), path(identities)
    path(auditor)
    output:
    tuple val(taxid), path('variants/variant_report.md'), emit: report
    tuple val(taxid), path('variants/*'), emit: files
    path('variants/versions.tsv'), emit: versions
    script:
    """
    set -euo pipefail
    mkdir variants
    bcftools query -l ${vcf} > variants/vcf_samples.txt
    bcftools index --stats ${vcf} > variants/coordinate_records.tsv
    bcftools stats ${vcf} > variants/bcftools_stats.txt
    python3 ${auditor} --directory variants --identities ${identities}
    printf 'process\\ttool\\tversion\\n' > variants/versions.tsv
    printf '%s\\tbcftools\\t%s\\n' '${task.process}' "\$(bcftools --version-only)" >> variants/versions.tsv
    """
}
