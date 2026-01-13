/*
========================================================================================
    FASTQC HiFi MODULE
========================================================================================
    Runs FastQC on HiFi reads
========================================================================================
*/

process FASTQC_HIFI {
    tag "${sample_id}"
    label 'fastqc'
    
    //temporarily output
    publishDir "${params.outdir}/${sample_id}/qc/hifi/fastqc", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hifi_fastq)
    
    output:
    tuple val(sample_id), path("*.html"), emit: fastqc_html
    tuple val(sample_id), path("*.zip"),  emit: fastqc_zip
    
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