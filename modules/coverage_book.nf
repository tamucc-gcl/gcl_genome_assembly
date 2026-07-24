process COVERAGE_BOOK {
    tag "$meta.id"
    label 'coverage_book'

    publishDir "${params.outdir}/bam/hifi/final", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly), path(hifi_reads)
    path(coverage_book_script)
    
    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    tuple val(meta), path("*.cov.${params.coverage_bin_bp}.bw"),                           emit: bigwig
    tuple val(meta), path("*.coverage_book.pdf"),                    emit: pdf
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: meta.id
    //def bin_size = task.ext.bin_size ?: 1000
    //def min_len = task.ext.min_len ?: 1000000
    //def min_mapq = task.ext.min_mapq ?: 5
    """
    # Index assembly
    samtools faidx ${assembly}
    
    # Map HiFi reads and sort
    minimap2 \\
        -t ${task.cpus} \\
        -ax map-hifi \\
        ${assembly} \\
        ${hifi_reads} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -o ${prefix}.sorted.bam
    
    # Index BAM
    samtools index -@ ${task.cpus} ${prefix}.sorted.bam
    
    # Generate coverage BigWig
    bamCoverage \\
        -b ${prefix}.sorted.bam \\
        -o ${prefix}.cov.${params.coverage_bin_bp}.bw \\
        --binSize ${params.coverage_bin_bp} \\
        --numberOfProcessors ${task.cpus} \\
        --ignoreDuplicates \\
        --minMappingQuality ${params.coverage_min_mapq}
    
    # Generate coverage book PDF
    python3 ${coverage_book_script} \\
        --bw ${prefix}.cov.${params.coverage_bin_bp}.bw \\
        --fai ${assembly}.fai \\
        --out_pdf ${prefix}.coverage_book.pdf \\
        --bin_size ${params.coverage_bin_bp} \\
        --min_len ${params.coverage_min_bp}
    """
    
    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai
    touch ${prefix}.cov.${params.coverage_bin_bp}.bw
    touch ${prefix}.coverage_book.pdf
    """
}