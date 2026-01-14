/*
========================================================================================
    HI-C CONTACT MAP MODULE
========================================================================================
    Generates Hi-C contact maps from mapped BAM files using cooler/HiCExplorer
    Produces contact matrices at multiple resolutions
    Can process either raw or filtered BAMs
========================================================================================
*/

process HIC_CONTACT_MAP {
    tag "${haplotype_id}_${qc_label}"
    label 'hic_contact_map'
    
    publishDir "${params.outdir}/qc/hic_mapping/${qc_label}/${haplotype_id}/contact_maps", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(bam), path(bai), path(assembly_fasta), val(qc_label)
    
    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.cool"), emit: cool
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.mcool"), emit: mcool
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_*.png"), emit: contact_maps
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_contact_stats.txt"), emit: stats
    
    script:
    def resolutions = params.hic_resolutions ?: "1000000,500000,100000,50000,10000"
    def filter_params = qc_label == "raw" ? 
        "--min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30" : 
        "--min-mapq 1 --walks-policy all"
    def dedup_step = qc_label == "raw" ?
        "| pairtools dedup --mark-dups --output-stats ${haplotype_id}_${qc_label}_dedup_stats.txt" :
        ""
    """
    # Create chromosome sizes file
    samtools faidx ${assembly_fasta}
    cut -f1,2 ${assembly_fasta}.fai > chrom.sizes
    
    # Convert BAM to pairs format
    # For raw BAMs: filter and deduplicate
    # For filtered BAMs: minimal processing (already clean)
    pairtools parse \\
        ${filter_params} \\
        --chroms-path chrom.sizes \\
        --output-stats ${haplotype_id}_${qc_label}_parse_stats.txt \\
        ${bam} \\
        | pairtools sort \\
            --nproc ${task.cpus} \\
            --tmpdir . \\
        ${dedup_step} \\
        --output ${haplotype_id}_${qc_label}.pairs.gz
    
    # Create cooler file at base resolution
    cooler cload pairs \\
        -c1 2 -p1 3 -c2 4 -p2 5 \\
        chrom.sizes:10000 \\
        ${haplotype_id}_${qc_label}.pairs.gz \\
        ${haplotype_id}_${qc_label}.cool
    
    # Create multi-resolution cooler file
    cooler zoomify \\
        --balance \\
        --resolutions ${resolutions} \\
        --out ${haplotype_id}_${qc_label}.mcool \\
        ${haplotype_id}_${qc_label}.cool
    
    # Generate contact map plots at different resolutions
    for res in 1000000 500000 100000; do
        hicPlotMatrix \\
            --matrix ${haplotype_id}_${qc_label}.mcool::/resolutions/\${res} \\
            --outFileName ${haplotype_id}_${qc_label}_\${res}bp_contact_map.png \\
            --title "${haplotype_id} (${qc_label}) Hi-C Contact Map (\${res}bp)" \\
            --log1p \\
            --dpi 300 \\
            --colorMap RdYlBu_r
    done
    
    # Generate contact statistics
    cat > ${haplotype_id}_${qc_label}_contact_stats.txt <<EOF
# Hi-C Contact Map Statistics for ${haplotype_id} (${qc_label})
# Generated: \$(date)

=== BAM Type ===
${qc_label}

=== Parsing Statistics ===
\$(cat ${haplotype_id}_${qc_label}_parse_stats.txt)
EOF

    # Add deduplication stats only for raw BAMs
    if [ "${qc_label}" = "raw" ]; then
        cat >> ${haplotype_id}_${qc_label}_contact_stats.txt <<EOF

=== Deduplication Statistics ===
\$(cat ${haplotype_id}_${qc_label}_dedup_stats.txt)
EOF
    fi

    cat >> ${haplotype_id}_${qc_label}_contact_stats.txt <<EOF

=== Cooler Information ===
\$(cooler info ${haplotype_id}_${qc_label}.cool)

=== Multi-resolution Cooler ===
\$(cooler tree ${haplotype_id}_${qc_label}.mcool)
EOF
    """
    
    stub:
    """
    touch ${haplotype_id}_${qc_label}.cool
    touch ${haplotype_id}_${qc_label}.mcool
    touch ${haplotype_id}_${qc_label}_1000000bp_contact_map.png
    touch ${haplotype_id}_${qc_label}_500000bp_contact_map.png
    touch ${haplotype_id}_${qc_label}_100000bp_contact_map.png
    touch ${haplotype_id}_${qc_label}_contact_stats.txt
    """
}