/*
========================================================================================
    TRIM HI-C MODULE
========================================================================================
    Trims Hi-C paired-end reads using fastp with adapter detection and quality filtering
    Repo location: modules/trim_hic.nf
========================================================================================
*/

process TRIM_HIC {
    tag "${meta.sample}"
    label 'fastp'

    //temporarily publish
    publishDir "${params.outdir}/fastq/hic/trimmed", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hic_r1), path(hic_r2)

    output:
    tuple val(meta), path("${meta.sample}_R1.trim.fastq.gz"), path("${meta.sample}_R2.trim.fastq.gz"), emit: trimmed_reads
    tuple val(meta), path("${meta.sample}_fastp.html"), emit: html
    tuple val(meta), path("${meta.sample}_fastp.json"), emit: json

    script:
    """
    fastp \\
        -i ${hic_r1} -I ${hic_r2} \\
        -o ${meta.sample}_R1.trim.fastq.gz -O ${meta.sample}_R2.trim.fastq.gz \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --cut_tail \\
        --cut_tail_window_size 4 \\
        --cut_tail_mean_quality 20 \\
        --length_required 30 \\
        --thread ${task.cpus} \\
        --html ${meta.sample}_fastp.html \\
        --json ${meta.sample}_fastp.json
    """

    stub:
    """
    touch ${meta.sample}_R1.trim.fastq.gz
    touch ${meta.sample}_R2.trim.fastq.gz
    touch ${meta.sample}_fastp.html
    touch ${meta.sample}_fastp.json
    """
}
