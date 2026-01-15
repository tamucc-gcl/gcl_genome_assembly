/*
========================================================================================
    HI-C SCAFFOLDING MODULE
========================================================================================
    Scaffolds assembled contigs using Hi-C contact information
    Uses YaHS (Yet another Hi-C scaffolding tool) for efficient scaffolding
    
    Input: Filtered Hi-C BAM file (valid pairs only) and assembly FASTA
    Output: Scaffolded assembly, AGP file, and bin files for visualization
========================================================================================
*/

process SCAFFOLD_HIC {
    tag "${haplotype_id}"
    label 'scaffold_hic'
    
    publishDir "${params.outdir}/scaffolding/yahs/${haplotype_id}", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(bam), path(bai), path(assembly_fasta)
    
    output:
    tuple val(haplotype_id), path("${haplotype_id}_scaffolds.fa"), emit: scaffolds
    tuple val(haplotype_id), path("${haplotype_id}_scaffolds_final.agp"), emit: agp
    tuple val(haplotype_id), path("${haplotype_id}_scaffolds_final.bin"), emit: bin
    tuple val(haplotype_id), path("${haplotype_id}.yahs.log"), emit: log
    
    script:
    def min_contig_len = params.yahs_min_contig_length ?: 10000
    def min_mapq = params.yahs_min_mapq ?: 1
    def resolution = params.yahs_resolution ?: 50000
    """
    # Index assembly if needed
    if [[ ! -s "${assembly_fasta}.fai" ]]; then
        samtools faidx ${assembly_fasta}
    fi
    
    # Run YaHS scaffolding
    yahs \\
        ${assembly_fasta} \\
        ${bam} \\
        -o ${haplotype_id}_scaffolds \\
        -l ${min_contig_len} \\
        -q ${min_mapq} \\
        -r ${resolution} \\
        2>&1 | tee ${haplotype_id}.yahs.log
    
    # Rename output files to standardized names
    mv ${haplotype_id}_scaffolds_scaffolds_final.fa ${haplotype_id}_scaffolds.fa
    mv ${haplotype_id}_scaffolds_scaffolds_final.agp ${haplotype_id}_scaffolds_final.agp
    """
    
    stub:
    """
    touch ${haplotype_id}_scaffolds.fa
    touch ${haplotype_id}_scaffolds_final.agp
    touch ${haplotype_id}_scaffolds_final.bin
    touch ${haplotype_id}.yahs.log
    """
}