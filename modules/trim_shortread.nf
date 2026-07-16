/*
========================================================================================
    TRIM SHORT-READ MODULE
========================================================================================
    Repo location: modules/trim_shortread.nf

    Adapter + quality trimming of raw short-read (Illumina) paired-end shotgun reads with
    fastp. Mirrors TRIM_HIC (same 'fastp' label) — the short-read analogue of the Hi-C
    raw->trim step, so assembly can start from raw reads rather than an external
    pre-cleaned product.

    Deliberately NO deduplication: mito/chloroplast are high-copy and a WGS library is
    heavily over-covered for them, so read-level dedup would strip genuine organelle depth.
    Nuclear assembly tolerates PCR duplicates; contamination is handled downstream by
    FCS-GX / FCS-adaptor at the assembly level. (If a nuclear/genome-size-only dedup stream
    is ever wanted, add --dedup via params.shortread_fastp_extra.)
========================================================================================
*/

process TRIM_SHORTREAD {
    tag "${meta.sample}"
    label 'fastp'

    publishDir "${params.outdir}/fastq/shortread/trimmed", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sr_r1), path(sr_r2)

    output:
    tuple val(meta), path("${meta.sample}_R1.trim.fastq.gz"), path("${meta.sample}_R2.trim.fastq.gz"), emit: trimmed_reads
    tuple val(meta), path("${meta.sample}_fastp.html"), emit: html
    tuple val(meta), path("${meta.sample}_fastp.json"), emit: json

    script:
    def qual   = params.shortread_fastp_cut_tail_quality ?: 20
    def window = params.shortread_fastp_cut_tail_window   ?: 4
    def minlen = params.shortread_fastp_length_required   ?: 30
    def extra  = params.shortread_fastp_extra ?: ''
    """
    fastp \\
        -i ${sr_r1} -I ${sr_r2} \\
        -o ${meta.sample}_R1.trim.fastq.gz -O ${meta.sample}_R2.trim.fastq.gz \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --cut_tail \\
        --cut_tail_window_size ${window} \\
        --cut_tail_mean_quality ${qual} \\
        --length_required ${minlen} \\
        --thread ${task.cpus} \\
        --html ${meta.sample}_fastp.html \\
        --json ${meta.sample}_fastp.json ${extra}
    """

    stub:
    """
    touch ${meta.sample}_R1.trim.fastq.gz
    touch ${meta.sample}_R2.trim.fastq.gz
    touch ${meta.sample}_fastp.html
    touch ${meta.sample}_fastp.json
    """
}
