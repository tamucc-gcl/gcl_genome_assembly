/*
========================================================================================
    MAP HI-C TO ASSEMBLY MODULE
========================================================================================
    Maps Hi-C paired-end reads to each haplotype assembly
    Uses BWA MEM for alignment, followed by sorting and indexing
========================================================================================
*/

process MAP_HIC_TO_ASSEMBLY {
    tag "${haplotype_id}"
    label 'map_hic'
    
    publishDir "${params.outdir}/mapping/hic/raw", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly_fasta), path(hic_r1), path(hic_r2)
    
    output:
    tuple val(haplotype_id), path("${haplotype_id}.sorted.bam"), path("${haplotype_id}.sorted.bam.bai"), emit: bam
    tuple val(haplotype_id), path("${haplotype_id}_mapping_stats.txt"), emit: stats
    
    script:
    """
    # Index the assembly
    bwa index ${assembly_fasta}
    samtools faidx ${assembly_fasta}
    
    # Map Hi-C reads to assembly
    bwa mem \\
        -t ${task.cpus} \\
        -5SP \\
        ${assembly_fasta} \\
        ${hic_r1} \\
        ${hic_r2} \\
        | samtools view -@ ${task.cpus} -Sb - \\
        | samtools sort -@ ${task.cpus} -o ${haplotype_id}.sorted.bam -
    
    # Index BAM file
    samtools index -@ ${task.cpus} ${haplotype_id}.sorted.bam
    
    # Generate mapping statistics
    samtools flagstat ${haplotype_id}.sorted.bam > ${haplotype_id}_mapping_stats.txt
    """
    
    stub:
    """
    touch ${haplotype_id}.sorted.bam
    touch ${haplotype_id}.sorted.bam.bai
    touch ${haplotype_id}_mapping_stats.txt
    """
}