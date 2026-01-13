/*
========================================================================================
    TRIM HI-C MODULE
========================================================================================
    Trims Hi-C paired-end reads using fastp with adapter detection and quality filtering
========================================================================================
*/

process TRIM_HIC {
    tag "${sample_id}"
    label 'fastp'

    //temporarily publish
    publishDir "${params.outdir}/${sample_id}/trimmed/hic", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hic_r1), path(hic_r2)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.html"), emit: html
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    
    script:
    """
    fastp \\
        -i ${hic_r1} -I ${hic_r2} \\
        -o ${sample_id}_R1.trim.fastq.gz -O ${sample_id}_R2.trim.fastq.gz \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --cut_tail \\
        --cut_tail_window_size 4 \\
        --cut_tail_mean_quality 20 \\
        --length_required 30 \\
        --thread ${task.cpus} \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json
    """
    
    stub:
    """
    touch ${sample_id}_R1.trim.fastq.gz
    touch ${sample_id}_R2.trim.fastq.gz
    touch ${sample_id}_fastp.html
    touch ${sample_id}_fastp.json
    """
}