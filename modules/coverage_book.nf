process COVERAGE_BOOK {
    tag "$meta.id"
    label 'coverage_book'

    publishDir "${params.outdir}/qc/assembly/coverage", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(fai)
    path(coverage_book_script)
    
    output:
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
    set -euo pipefail
    
    # Generate coverage BigWig
    bamCoverage \\
        -b ${bam} \\
        -o ${prefix}.cov.${params.coverage_bin_bp}.bw \\
        --binSize ${params.coverage_bin_bp} \\
        --numberOfProcessors ${task.cpus} \\
        --ignoreDuplicates \\
        --minMappingQuality ${params.coverage_min_mapq}
    
    # Generate coverage book PDF
    python3 ${coverage_book_script} \\
        --bw ${prefix}.cov.${params.coverage_bin_bp}.bw \\
        --fai ${fai} \\
        --out_pdf ${prefix}.coverage_book.pdf \\
        --bin_size ${params.coverage_bin_bp} \\
        --min_len ${params.coverage_min_bp}
    """
    
    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.cov.${params.coverage_bin_bp}.bw
    touch ${prefix}.coverage_book.pdf
    """
}
