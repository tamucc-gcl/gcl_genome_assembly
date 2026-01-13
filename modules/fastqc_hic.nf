/*
========================================================================================
    FASTQC HI-C MODULE
========================================================================================
    Runs FastQC on Hi-C paired-end reads
========================================================================================
*/

process FASTQC_HIC {
    tag "${sample_id}"
    label 'fastqc'
    
    //temporarily output
    publishDir "${params.outdir}/${sample_id}/qc/hic/${qc_label}/fastqc", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hic_r1), path(hic_r2)
    val qc_label  // "raw" or "trimmed"
    
    output:
    tuple val(sample_id), path("*.html"), emit: fastqc_html
    tuple val(sample_id), path("*.zip"),  emit: fastqc_zip
    
    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${hic_r1} \\
        ${hic_r2}
    """
    
    stub:
    """
    touch ${hic_r1.baseName}_fastqc.html
    touch ${hic_r1.baseName}_fastqc.zip
    touch ${hic_r2.baseName}_fastqc.html
    touch ${hic_r2.baseName}_fastqc.zip
    """
}