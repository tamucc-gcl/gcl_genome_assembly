/*
========================================================================================
    BAM TO FASTQ MODULE
========================================================================================
    Converts HiFi BAM files to compressed FASTQ format using samtools and pigz
========================================================================================
*/

process BAM_TO_FASTQ {
    tag "${sample_id}"
    label 'bam_to_fastq'
    
    //temporarily publish for debugging
    publishDir "${params.outdir}/${sample_id}/fastq", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hifi_bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz"), emit: fastq
    
    script:
    """
    samtools fastq \\
        -@ ${task.cpus} \\
        ${hifi_bam} \\
        | pigz -p ${task.cpus} > ${sample_id}.fastq.gz
    """
    
    stub:
    """
    touch ${sample_id}.fastq.gz
    """
}