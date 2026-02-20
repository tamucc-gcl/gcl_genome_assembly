/*
========================================================================================
    PURGE_DUPS MODULE
========================================================================================
    Remove haplotig duplications from genome assemblies using purge_dups
    
    This module performs:
    1. Map HiFi reads to assembly (minimap2)
    2. Calculate coverage statistics (pbcstat)
    3. Determine coverage cutoffs (calcuts)
    4. Self-alignment for duplicate detection (minimap2)
    5. Identify duplications (purge_dups)
    6. Generate purged assembly and haplotigs (get_seqs)
    
    Inputs:
    - haplotype_id: Unique identifier for the haplotype (e.g., sample_hap1)
    - assembly: FASTA file of the assembly to purge
    - hifi_reads: HiFi reads FASTQ file for coverage calculation
    
    Outputs:
    - purged_assembly: Assembly with haplotigs removed
    - haplotigs: Removed haplotig sequences
    - coverage_stats: Coverage statistics and cutoffs for QC
========================================================================================
*/

process PURGE_DUPS {
    tag "${haplotype_id}"
    label 'purge_dups'
    
    publishDir "${params.outdir}/contig/purge_dups", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly), path(hifi_reads)
    
    output:
    tuple val(haplotype_id), path("${haplotype_id}.purged.fa"), emit: purged_assembly
    tuple val(haplotype_id), path("${haplotype_id}.haplotigs.fa"), emit: haplotigs
    tuple val(haplotype_id), path("${haplotype_id}.purge_dups_stats"), emit: stats_dir
    tuple val(haplotype_id), path("${haplotype_id}.purge_dups.log"), emit: log
    
    script:
    def minimap_threads = Math.max(1, task.cpus - 4)
    def pigz_threads = Math.min(4, task.cpus)
    """
    set -euo pipefail
    
    # Create stats output directory
    mkdir -p ${haplotype_id}.purge_dups_stats
    
    echo "=== Step 1: Mapping HiFi reads to assembly ===" | tee ${haplotype_id}.purge_dups.log
    minimap2 -xmap-hifi -t ${minimap_threads} ${assembly} ${hifi_reads} | \\
        pigz -p ${pigz_threads} -c > aligned.paf.gz
    
    echo "=== Step 2: Calculating coverage statistics ===" | tee -a ${haplotype_id}.purge_dups.log
    pbcstat aligned.paf.gz 2>&1 | tee -a ${haplotype_id}.purge_dups.log
    
    # Copy coverage stats to output directory
    cp PB.stat ${haplotype_id}.purge_dups_stats/
    cp PB.base.cov ${haplotype_id}.purge_dups_stats/
    
    echo "=== Step 3: Calculating coverage cutoffs ===" | tee -a ${haplotype_id}.purge_dups.log
    calcuts PB.stat > cutoffs 2> calcuts.log
    cat calcuts.log | tee -a ${haplotype_id}.purge_dups.log
    cp cutoffs ${haplotype_id}.purge_dups_stats/
    cp calcuts.log ${haplotype_id}.purge_dups_stats/
    
    echo "=== Step 4: Splitting assembly for self-alignment ===" | tee -a ${haplotype_id}.purge_dups.log
    split_fa ${assembly} > split.fa
    
    echo "=== Step 5: Self-alignment ===" | tee -a ${haplotype_id}.purge_dups.log
    minimap2 -xasm5 -DP -t ${minimap_threads} split.fa split.fa | \\
        pigz -p ${pigz_threads} -c > split.self.paf.gz
    
    echo "=== Step 6: Identifying duplications ===" | tee -a ${haplotype_id}.purge_dups.log
    purge_dups -2 -T cutoffs -c PB.base.cov split.self.paf.gz > dups.bed 2>&1 | \\
        tee -a ${haplotype_id}.purge_dups.log
    cp dups.bed ${haplotype_id}.purge_dups_stats/
    
    echo "=== Step 7: Generating purged assembly ===" | tee -a ${haplotype_id}.purge_dups.log
    get_seqs -e dups.bed ${assembly}
    
    # Rename outputs to include haplotype_id
    mv purged.fa ${haplotype_id}.purged.fa
    mv hap.fa ${haplotype_id}.haplotigs.fa
    
    # Calculate and log summary statistics
    echo "" | tee -a ${haplotype_id}.purge_dups.log
    echo "=== Summary ===" | tee -a ${haplotype_id}.purge_dups.log
    
    # Original assembly size
    orig_size=\$(awk '/^>/{if(s)print s; s=0; next} {s+=length}END{print s}' ${assembly})
    echo "Original assembly size: \${orig_size}" | tee -a ${haplotype_id}.purge_dups.log
    
    # Purged assembly size
    purged_size=\$(awk '/^>/{if(s)print s; s=0; next} {s+=length}END{print s}' ${haplotype_id}.purged.fa)
    echo "Purged assembly size: \${purged_size}" | tee -a ${haplotype_id}.purge_dups.log
    
    # Haplotigs size
    hap_size=\$(awk '/^>/{if(s)print s; s=0; next} {s+=length}END{print s}' ${haplotype_id}.haplotigs.fa)
    echo "Haplotigs removed: \${hap_size}" | tee -a ${haplotype_id}.purge_dups.log
    
    # Contigs counts
    orig_n=\$(grep -c '^>' ${assembly})
    purged_n=\$(grep -c '^>' ${haplotype_id}.purged.fa)
    hap_n=\$(grep -c '^>' ${haplotype_id}.haplotigs.fa)
    echo "Original contigs: \${orig_n}" | tee -a ${haplotype_id}.purge_dups.log
    echo "Purged contigs: \${purged_n}" | tee -a ${haplotype_id}.purge_dups.log
    echo "Haplotig contigs: \${hap_n}" | tee -a ${haplotype_id}.purge_dups.log
    
    # Save summary to stats directory
    cat > ${haplotype_id}.purge_dups_stats/summary.tsv <<EOF
haplotype_id\toriginal_size\tpurged_size\thaplotigs_size\toriginal_contigs\tpurged_contigs\thaplotig_contigs
${haplotype_id}\t\${orig_size}\t\${purged_size}\t\${hap_size}\t\${orig_n}\t\${purged_n}\t\${hap_n}
EOF
    
    echo "=== purge_dups completed successfully ===" | tee -a ${haplotype_id}.purge_dups.log
    
    # Cleanup intermediate files
    rm -f aligned.paf.gz split.fa split.self.paf.gz PB.stat PB.base.cov cutoffs dups.bed
    """
    
    stub:
    """
    mkdir -p ${haplotype_id}.purge_dups_stats
    touch ${haplotype_id}.purged.fa
    touch ${haplotype_id}.haplotigs.fa
    touch ${haplotype_id}.purge_dups.log
    touch ${haplotype_id}.purge_dups_stats/summary.tsv
    touch ${haplotype_id}.purge_dups_stats/PB.stat
    touch ${haplotype_id}.purge_dups_stats/cutoffs
    touch ${haplotype_id}.purge_dups_stats/dups.bed
    """
}