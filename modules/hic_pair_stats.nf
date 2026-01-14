/*
========================================================================================
    HI-C PAIR STATISTICS MODULE
========================================================================================
    Analyzes Hi-C mapping pairs for quality metrics:
    - Valid pairs percentage
    - Trans/cis ratios
    - Insert size distribution
    - Mapping quality distribution
    Can process either raw or filtered BAMs
========================================================================================
*/

process HIC_PAIR_STATS {
    tag "${haplotype_id}_${qc_label}"
    label 'hic_pair_stats'
    
    publishDir "${params.outdir}/qc/hic_mapping/${qc_label}/${haplotype_id}/pair_stats", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(bam), path(bai), val(qc_label)
    
    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_pair_types.txt"), emit: pair_types
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_trans_cis_ratio.txt"), emit: trans_cis
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_insert_size_dist.txt"), emit: insert_dist
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_pair_stats_summary.txt"), emit: summary
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_*.pdf"), emit: plots
    
    script:
    """
    # Extract pair statistics using samtools
    # Get mapping quality distribution
    samtools view -F 256 ${bam} \\
        | awk '{print \$5}' \\
        | sort -n \\
        | uniq -c \\
        > ${haplotype_id}_${qc_label}_mapq_dist.txt
    
    # Analyze pair types and orientations
    samtools view -f 1 -F 256 ${bam} \\
        | awk '{
            # Get flags
            flag = \$2
            # Determine orientation
            r1_rev = and(flag, 16) ? "-" : "+"
            r2_rev = and(flag, 32) ? "-" : "+"
            
            # Same chromosome (cis) or different (trans)
            if (\$3 == \$7) {
                pair_type = "cis"
                # Calculate insert size for cis pairs
                if (\$9 != 0) {
                    insert = (\$9 < 0) ? -\$9 : \$9
                    if (insert < 1000) bin = "0-1kb"
                    else if (insert < 10000) bin = "1-10kb"
                    else if (insert < 100000) bin = "10-100kb"
                    else if (insert < 1000000) bin = "100kb-1Mb"
                    else bin = ">1Mb"
                    print "cis_distance", bin
                }
            } else {
                pair_type = "trans"
            }
            
            # Count orientations
            orientation = r1_rev r2_rev
            print pair_type, orientation
        }' \\
        | sort \\
        | uniq -c \\
        > ${haplotype_id}_${qc_label}_pair_types.txt
    
    # Calculate trans/cis ratio
    awk '
        /cis/ {cis += \$1}
        /trans/ {trans += \$1}
        END {
            total = cis + trans
            if (total > 0) {
                print "Total pairs:", total
                print "Cis pairs:", cis, sprintf("(%.2f%%)", 100*cis/total)
                print "Trans pairs:", trans, sprintf("(%.2f%%)", 100*trans/total)
                if (cis > 0) print "Trans/Cis ratio:", trans/cis
            }
        }
    ' ${haplotype_id}_${qc_label}_pair_types.txt > ${haplotype_id}_${qc_label}_trans_cis_ratio.txt
    
    # Get insert size distribution for cis pairs
    samtools view -f 1 -F 256 ${bam} \\
        | awk '\$3 == \$7 && \$9 > 0 {print \$9}' \\
        | sort -n \\
        | uniq -c \\
        > ${haplotype_id}_${qc_label}_insert_size_dist.txt
    
    # Create comprehensive summary
    cat > ${haplotype_id}_${qc_label}_pair_stats_summary.txt <<EOF
# Hi-C Pair Statistics Summary for ${haplotype_id} (${qc_label})
# Generated: \$(date)

=== BAM Type ===
${qc_label}

=== Basic Mapping Statistics ===
\$(samtools flagstat ${bam})

=== Pair Type Distribution ===
\$(cat ${haplotype_id}_${qc_label}_pair_types.txt)

=== Trans/Cis Analysis ===
\$(cat ${haplotype_id}_${qc_label}_trans_cis_ratio.txt)

=== Insert Size Statistics (cis pairs) ===
Mean: \$(awk '{sum+=\$2*\$1; n+=\$1} END {if(n>0) print sum/n; else print 0}' ${haplotype_id}_${qc_label}_insert_size_dist.txt)
Median: \$(awk '{for(i=0;i<\$1;i++) print \$2}' ${haplotype_id}_${qc_label}_insert_size_dist.txt | sort -n | awk '{a[NR]=\$0} END {print a[int(NR/2)]}')

=== Mapping Quality Distribution ===
\$(head -20 ${haplotype_id}_${qc_label}_mapq_dist.txt)
EOF

    # Generate plots using R
    Rscript - <<'RSCRIPT'
    library(ggplot2)
    library(gridExtra)
    
    # Read data
    pair_types <- read.table("${haplotype_id}_${qc_label}_pair_types.txt")
    colnames(pair_types) <- c("count", "type", "orientation")
    
    insert_dist <- read.table("${haplotype_id}_${qc_label}_insert_size_dist.txt")
    colnames(insert_dist) <- c("count", "insert_size")
    
    mapq_dist <- read.table("${haplotype_id}_${qc_label}_mapq_dist.txt")
    colnames(mapq_dist) <- c("count", "mapq")
    
    # Plot 1: Pair type distribution
    p1 <- ggplot(pair_types, aes(x=type, y=count, fill=orientation)) +
        geom_bar(stat="identity") +
        theme_minimal() +
        labs(title="Hi-C Pair Types (${qc_label})", x="Pair Type", y="Count") +
        theme(axis.text.x = element_text(angle=45, hjust=1))
    
    # Plot 2: Insert size distribution (log scale)
    insert_filtered <- insert_dist[insert_dist\$insert_size < 1000000,]
    p2 <- ggplot(insert_filtered, aes(x=insert_size, y=count)) +
        geom_line() +
        scale_x_log10() +
        theme_minimal() +
        labs(title="Insert Size Distribution (${qc_label})", 
             x="Insert Size (bp, log scale)", y="Count")
    
    # Plot 3: Mapping quality distribution
    p3 <- ggplot(mapq_dist, aes(x=mapq, y=count)) +
        geom_bar(stat="identity", fill="steelblue") +
        theme_minimal() +
        labs(title="Mapping Quality Distribution (${qc_label})", x="MAPQ", y="Count")
    
    # Save combined plot
    pdf("${haplotype_id}_${qc_label}_hic_pair_stats.pdf", width=12, height=10)
    grid.arrange(p1, p2, p3, ncol=2)
    dev.off()
RSCRIPT
    """
    
    stub:
    """
    touch ${haplotype_id}_${qc_label}_pair_types.txt
    touch ${haplotype_id}_${qc_label}_trans_cis_ratio.txt
    touch ${haplotype_id}_${qc_label}_insert_size_dist.txt
    touch ${haplotype_id}_${qc_label}_pair_stats_summary.txt
    touch ${haplotype_id}_${qc_label}_hic_pair_stats.pdf
    """
}