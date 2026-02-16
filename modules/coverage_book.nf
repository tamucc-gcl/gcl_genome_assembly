process COVERAGE_BOOK {
    tag "$meta.id"
    label 'coverage_book'

    publishDir "${params.outdir}/bam/hifi/final", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly), path(hifi_reads)
    
    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    tuple val(meta), path("*.cov.1kb.bw"),                           emit: bigwig
    tuple val(meta), path("*.coverage_book.pdf"),                    emit: pdf
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: meta.id
    //def bin_size = task.ext.bin_size ?: 1000
    //def min_len = task.ext.min_len ?: 1000000
    //def min_mapq = task.ext.min_mapq ?: 5
    """
    bin_size=1000
    min_len=1000000
    min_mapq=5

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
        -o ${prefix}.cov.1kb.bw \\
        --binSize ${bin_size} \\
        --numberOfProcessors ${task.cpus} \\
        --ignoreDuplicates \\
        --minMappingQuality ${min_mapq}
    
    # Generate coverage book PDF
    python3 ${projectDir}/py_scripts/bigwig_genome_book.py \\
        --bw ${prefix}.cov.1kb.bw \\
        --fai ${assembly}.fai \\
        --out_pdf ${prefix}.coverage_book.pdf \\
        --bin_size ${bin_size} \\
        --min_len ${min_len}
    """
    
    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai
    touch ${prefix}.cov.1kb.bw
    touch ${prefix}.coverage_book.pdf
    """
}