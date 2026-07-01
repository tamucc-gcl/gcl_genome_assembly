/*
========================================================================================
    BAM TO FASTQ MODULE
========================================================================================
    Converts HiFi BAM files to compressed FASTQ format using samtools and pigz
    Repo location: modules/bam_to_fastq.nf
========================================================================================
*/

process BAM_TO_FASTQ {
    tag "${meta.sample}"
    label 'bam_to_fastq'

    //temporarily publish for debugging
    publishDir "${params.outdir}/fastq/hifi", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hifi_bam)

    output:
    tuple val(meta), path("${meta.sample}.fastq.gz"), emit: fastq

    script:
    """
    samtools fastq \\
        -@ ${task.cpus} \\
        ${hifi_bam} \\
        | pigz -p ${task.cpus} > ${meta.sample}.fastq.gz
    """

    stub:
    """
    touch ${meta.sample}.fastq.gz
    """
}
