/*
========================================================================================
    HI-C COVERAGE MODULE
========================================================================================
    Calculates Hi-C read coverage across the assembly
    Identifies regions with low/high coverage and potential issues
    Can process either raw or filtered BAMs
========================================================================================
*/

process HIC_COVERAGE {
    tag "${haplotype_id}_${qc_label}"
    label 'hic_coverage'
    
    publishDir "${params.outdir}/qc/hic_mapping/${qc_label}/${haplotype_id}/coverage", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(bam), path(bai), path(assembly_fasta), val(qc_label)
    
    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_coverage.txt"), emit: coverage
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_coverage_stats.txt"), emit: stats
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_coverage_*.png"), emit: plots
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.bedgraph"), emit: bedgraph
    
    script:
    def window_size = params.hic_coverage_window ?: 100000
    """
    # Calculate coverage in windows
    samtools faidx ${assembly_fasta}
    
    # Create windows across the genome
    bedtools makewindows \\
        -g ${assembly_fasta}.fai \\
        -w ${window_size} \\
        > windows.bed
    
    # Calculate coverage per window
    bedtools coverage \\
        -a windows.bed \\
        -b ${bam} \\
        -sorted \\
        > ${haplotype_id}_${qc_label}_coverage.txt
    
    # Generate bedgraph for visualization
    samtools depth -a ${bam} \\
        | awk '{print \$1"\\t"\$2-1"\\t"\$2"\\t"\$3}' \\
        > ${haplotype_id}_${qc_label}.bedgraph
    
    # Calculate coverage statistics
    awk '{
        count++
        sum += \$4
        if (count == 1 || \$4 < min) min = \$4
        if (count == 1 || \$4 > max) max = \$4
        cov[\$4]++
    } END {
        mean = sum / count
        print "BAM Type: ${qc_label}"
        print "Total windows:", count
        print "Mean coverage:", mean
        print "Min coverage:", min
        print "Max coverage:", max
        print "Coverage distribution:"
        for (c in cov) {
            print "  Coverage", c":", cov[c], "windows"
        }
    }' ${haplotype_id}_${qc_label}_coverage.txt > ${haplotype_id}_${qc_label}_coverage_stats.txt
    
    # Generate coverage plots using R
    Rscript - <<'RSCRIPT'
    library(ggplot2)
    library(dplyr)
    
    # Read coverage data
    cov <- read.table("${haplotype_id}_${qc_label}_coverage.txt", 
                      col.names=c("chr", "start", "end", "count", "bases", "length", "fraction"))
    
    # Calculate normalized coverage
    cov\$normalized_cov <- cov\$count / (cov\$length / 1000)  # reads per kb
    
    # Get chromosome sizes for proper ordering
    chr_sizes <- cov %>% 
        group_by(chr) %>% 
        summarize(max_pos = max(end)) %>%
        arrange(desc(max_pos))
    
    cov\$chr <- factor(cov\$chr, levels=chr_sizes\$chr)
    
    # Plot 1: Coverage across genome
    p1 <- ggplot(cov, aes(x=start/1e6, y=normalized_cov)) +
        geom_line(alpha=0.7) +
        facet_wrap(~chr, scales="free_x", ncol=2) +
        theme_minimal() +
        labs(title=paste("Hi-C Coverage (${qc_label}) -", "${haplotype_id}"),
             x="Position (Mb)", y="Normalized Coverage (reads/kb)") +
        theme(strip.text = element_text(size=8))
    
    ggsave("${haplotype_id}_${qc_label}_coverage_genome.png", p1, width=14, height=10, dpi=300)
    
    # Plot 2: Coverage distribution
    p2 <- ggplot(cov, aes(x=normalized_cov)) +
        geom_histogram(bins=50, fill="steelblue", alpha=0.7) +
        scale_x_log10() +
        theme_minimal() +
        labs(title=paste("Hi-C Coverage Distribution (${qc_label})"),
             x="Normalized Coverage (reads/kb, log scale)", y="Count")
    
    ggsave("${haplotype_id}_${qc_label}_coverage_distribution.png", p2, width=8, height=6, dpi=300)
    
    # Plot 3: Cumulative coverage
    cov_sorted <- sort(cov\$normalized_cov)
    cumulative <- data.frame(
        coverage = cov_sorted,
        percentile = (1:length(cov_sorted)) / length(cov_sorted) * 100
    )
    
    p3 <- ggplot(cumulative, aes(x=coverage, y=percentile)) +
        geom_line(color="darkblue", size=1) +
        theme_minimal() +
        labs(title=paste("Cumulative Hi-C Coverage (${qc_label})"),
             x="Normalized Coverage (reads/kb)", y="Percentile (%)") +
        geom_vline(xintercept=median(cov\$normalized_cov), 
                   linetype="dashed", color="red", alpha=0.7) +
        annotate("text", x=median(cov\$normalized_cov)*1.5, y=50, 
                 label=paste("Median:", round(median(cov\$normalized_cov), 2)))
    
    ggsave("${haplotype_id}_${qc_label}_coverage_cumulative.png", p3, width=8, height=6, dpi=300)
RSCRIPT
    """
    
    stub:
    """
    touch ${haplotype_id}_${qc_label}_coverage.txt
    touch ${haplotype_id}_${qc_label}_coverage_stats.txt
    touch ${haplotype_id}_${qc_label}_coverage_genome.png
    touch ${haplotype_id}_${qc_label}_coverage_distribution.png
    touch ${haplotype_id}_${qc_label}_coverage_cumulative.png
    touch ${haplotype_id}_${qc_label}.bedgraph
    """
}