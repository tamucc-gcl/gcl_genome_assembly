/*
========================================================================================
    FASTQC HI-C MODULE
========================================================================================
    Repo location: modules/fastqc_hic.nf
    Runs FastQC on Hi-C paired-end reads (per-sample; carries meta, meta.id == sample).
========================================================================================
*/

process FASTQC_HIC {
    tag "${meta.id}"
    label 'fastqc'

    //temporarily output
    //publishDir "${params.outdir}/${meta.id}/qc/hic/${qc_label}/fastqc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hic_r1), path(hic_r2)
    val qc_label  // "raw" or "trimmed"

    output:
    tuple val(meta), path("*.html"), emit: fastqc_html
    tuple val(meta), path("*.zip"),  emit: fastqc_zip
    path "versions.tsv", emit: versions

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${hic_r1} \\
        ${hic_r2}

    printf 'FastQC\t%s\n' "$(fastqc --version 2>&1 | sed 's/FastQC //')" > versions.tsv
    """

    stub:
    """
    touch ${hic_r1.baseName}_fastqc.html
    touch ${hic_r1.baseName}_fastqc.zip
    touch ${hic_r2.baseName}_fastqc.html
    touch ${hic_r2.baseName}_fastqc.zip
    touch versions.tsv
    """
}
