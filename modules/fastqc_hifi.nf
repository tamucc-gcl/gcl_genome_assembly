/*
========================================================================================
    FASTQC HiFi MODULE
========================================================================================
    Repo location: modules/fastqc_hifi.nf
    Runs FastQC on HiFi reads (per-sample; carries meta, meta.id == sample).
========================================================================================
*/

process FASTQC_HIFI {
    tag "${meta.id}"
    label 'fastqc'

    //temporarily output
    //publishDir "${params.outdir}/${meta.id}/qc/hifi/fastqc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hifi_fastq)

    output:
    tuple val(meta), path("*.html"), emit: fastqc_html
    tuple val(meta), path("*.zip"),  emit: fastqc_zip

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${hifi_fastq}
    """

    stub:
    """
    touch ${hifi_fastq.baseName}_fastqc.html
    touch ${hifi_fastq.baseName}_fastqc.zip
    """
}
