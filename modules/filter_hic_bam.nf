/*
========================================================================================
    FILTER HI-C BAM MODULE
========================================================================================
    Filters Hi-C BAM files to remove invalid pairs and PCR duplicates
    - Removes unmapped reads
    - Removes low MAPQ reads (< 30)
    - Removes improper pairs
    - Removes PCR duplicates using pairtools
    - Produces clean BAM files for downstream analysis
========================================================================================
*/

process FILTER_HIC_BAM {
    tag "${haplotype_id}"
    label 'filter_hic_bam'
    
    publishDir "${params.outdir}/mapping/hic/filtered", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(bam), path(bai), path(assembly_fasta)
    
    output:
    tuple val(haplotype_id), path("${haplotype_id}.filtered.sorted.bam"), path("${haplotype_id}.filtered.sorted.bam.bai"), emit: bam
    tuple val(haplotype_id), path("${haplotype_id}_filtering_stats.txt"), emit: stats
    tuple val(haplotype_id), path("${haplotype_id}.pairs.gz"), emit: pairs
    
    script:
    """
    # Create chromosome sizes file
    samtools faidx ${assembly_fasta}
    cut -f1,2 ${assembly_fasta}.fai > chrom.sizes
    
    # Convert BAM to pairs format, filtering and deduplicating in one pass
    pairtools parse \\
        --min-mapq ${params.hic_min_mapq ?: 30} \\
        --walks-policy 5unique \\
        --max-inter-align-gap 30 \\
        --chroms-path chrom.sizes \\
        --output-stats ${haplotype_id}_parse_stats.txt \\
        ${bam} \\
        | pairtools sort \\
            --nproc ${task.cpus} \\
            --tmpdir . \\
        | pairtools dedup \\
            --mark-dups \\
            --output-stats ${haplotype_id}_dedup_stats.txt \\
            --output-dups ${haplotype_id}.dups.pairs.gz \\
            --output ${haplotype_id}.pairs.gz
    
    # Convert filtered pairs back to BAM
    # Split pairs into two SAM files (R1 and R2)
    pairtools split \\
        --output-sam ${haplotype_id}.R1.sam ${haplotype_id}.R2.sam \\
        ${haplotype_id}.pairs.gz
    
    # Convert SAM to BAM, merge, sort and index
    samtools view -bS ${haplotype_id}.R1.sam > ${haplotype_id}.R1.bam
    samtools view -bS ${haplotype_id}.R2.sam > ${haplotype_id}.R2.bam
    
    # Merge R1 and R2 BAMs
    samtools merge -@ ${task.cpus} ${haplotype_id}.merged.bam ${haplotype_id}.R1.bam ${haplotype_id}.R2.bam
    
    # Sort and index
    samtools sort -@ ${task.cpus} -o ${haplotype_id}.filtered.sorted.bam ${haplotype_id}.merged.bam
    samtools index -@ ${task.cpus} ${haplotype_id}.filtered.sorted.bam
    
    # Generate comprehensive filtering statistics
    cat > ${haplotype_id}_filtering_stats.txt <<EOF
================================================================================
Hi-C BAM Filtering Statistics for ${haplotype_id}
================================================================================
Generated: \$(date)

================================================================================
ORIGINAL BAM STATISTICS
================================================================================
\$(samtools flagstat ${bam})

================================================================================
PARSING STATISTICS
================================================================================
\$(cat ${haplotype_id}_parse_stats.txt)

================================================================================
DEDUPLICATION STATISTICS
================================================================================
\$(cat ${haplotype_id}_dedup_stats.txt)

================================================================================
FILTERED BAM STATISTICS
================================================================================
\$(samtools flagstat ${haplotype_id}.filtered.sorted.bam)

================================================================================
FILTERING SUMMARY
================================================================================
EOF

    # Calculate filtering metrics
    original_reads=\$(samtools view -c -F 256 ${bam})
    filtered_reads=\$(samtools view -c -F 256 ${haplotype_id}.filtered.sorted.bam)
    retained_pct=\$(echo "scale=2; \$filtered_reads / \$original_reads * 100" | bc)
    
    cat >> ${haplotype_id}_filtering_stats.txt <<EOF
Original reads (excluding secondary): \${original_reads}
Filtered reads (valid pairs): \${filtered_reads}
Retention rate: \${retained_pct}%

Filtering criteria applied:
- Minimum MAPQ: ${params.hic_min_mapq ?: 30}
- Walks policy: 5unique (removes walks across multiple ligation junctions)
- Max inter-alignment gap: 30bp
- PCR duplicates: removed

Notes:
- Valid pairs are properly mapped Hi-C contacts
- Low retention (<50%) may indicate library quality issues
- High retention (>70%) typically indicates good Hi-C library
================================================================================
EOF

    # Clean up intermediate files
    rm -f ${haplotype_id}.R1.sam ${haplotype_id}.R2.sam
    rm -f ${haplotype_id}.R1.bam ${haplotype_id}.R2.bam
    rm -f ${haplotype_id}.merged.bam
    """
    
    stub:
    """
    touch ${haplotype_id}.filtered.sorted.bam
    touch ${haplotype_id}.filtered.sorted.bam.bai
    touch ${haplotype_id}_filtering_stats.txt
    touch ${haplotype_id}.pairs.gz
    """
}