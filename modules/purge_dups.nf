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
    - log: Processing log
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
    tuple val(haplotype_id), path("${haplotype_id}.purge_dups.log"), emit: log
    
    script:
    def minimap_threads = Math.max(1, task.cpus - 4)
    def pigz_threads = Math.min(4, task.cpus)
    """
    set -euo pipefail
    
    echo "=== purge_dups started ===" | tee ${haplotype_id}.purge_dups.log
    echo "Working directory: \$(pwd)" | tee -a ${haplotype_id}.purge_dups.log
    echo "" | tee -a ${haplotype_id}.purge_dups.log
    
    # Copy assembly to local file to avoid potential symlink issues with get_seqs
    echo "Copying assembly to local file..." | tee -a ${haplotype_id}.purge_dups.log
    cp ${assembly} assembly.fa
    
    echo "=== Step 1: Mapping HiFi reads to assembly ===" | tee -a ${haplotype_id}.purge_dups.log
    minimap2 -xmap-hifi -t ${minimap_threads} assembly.fa ${hifi_reads} | \\
        pigz -p ${pigz_threads} -c > aligned.paf.gz
    echo "Mapping complete" | tee -a ${haplotype_id}.purge_dups.log
    
    echo "=== Step 2: Calculating coverage statistics ===" | tee -a ${haplotype_id}.purge_dups.log
    pbcstat aligned.paf.gz 2>&1 | tee -a ${haplotype_id}.purge_dups.log
    
    echo "=== Step 3: Calculating coverage cutoffs ===" | tee -a ${haplotype_id}.purge_dups.log
    calcuts PB.stat > cutoffs 2>> ${haplotype_id}.purge_dups.log
    echo "Cutoffs:" | tee -a ${haplotype_id}.purge_dups.log
    cat cutoffs | tee -a ${haplotype_id}.purge_dups.log
    echo "" | tee -a ${haplotype_id}.purge_dups.log
    
    echo "=== Step 4: Splitting assembly for self-alignment ===" | tee -a ${haplotype_id}.purge_dups.log
    split_fa assembly.fa > split.fa
    echo "Split into \$(grep -c '^>' split.fa) sequences" | tee -a ${haplotype_id}.purge_dups.log
    
    echo "=== Step 5: Self-alignment ===" | tee -a ${haplotype_id}.purge_dups.log
    minimap2 -xasm5 -DP -t ${minimap_threads} split.fa split.fa | \\
        pigz -p ${pigz_threads} -c > split.self.paf.gz
    echo "Self-alignment complete" | tee -a ${haplotype_id}.purge_dups.log
    
    echo "=== Step 6: Identifying duplications ===" | tee -a ${haplotype_id}.purge_dups.log
    # CRITICAL: stdout goes to dups.bed, stderr (log messages) goes to log file
    # Do NOT use 2>&1 here or it corrupts the BED file!
    purge_dups -2 -T cutoffs -c PB.base.cov split.self.paf.gz > dups.bed 2>> ${haplotype_id}.purge_dups.log
    
    echo "dups.bed: \$(wc -l < dups.bed) lines" | tee -a ${haplotype_id}.purge_dups.log
    echo "First 10 lines of dups.bed:" | tee -a ${haplotype_id}.purge_dups.log
    head -10 dups.bed | tee -a ${haplotype_id}.purge_dups.log
    echo "" | tee -a ${haplotype_id}.purge_dups.log
    
    echo "=== Step 7: Generating purged assembly ===" | tee -a ${haplotype_id}.purge_dups.log
    
    if [ -s dups.bed ] && [ \$(wc -l < dups.bed) -gt 0 ]; then
        echo "Running: get_seqs -e dups.bed assembly.fa" | tee -a ${haplotype_id}.purge_dups.log
        get_seqs -e dups.bed assembly.fa 2>&1 | tee -a ${haplotype_id}.purge_dups.log
        echo "get_seqs completed" | tee -a ${haplotype_id}.purge_dups.log
    else
        echo "No duplications identified, using original assembly" | tee -a ${haplotype_id}.purge_dups.log
        cp assembly.fa purged.fa
        touch hap.fa
    fi
    
    # Calculate summary statistics
    echo "" | tee -a ${haplotype_id}.purge_dups.log
    echo "=== Summary ===" | tee -a ${haplotype_id}.purge_dups.log
    
    orig_size=\$(awk '/^>/{if(s)print s; s=0; next} {s+=length}END{print s}' assembly.fa)
    purged_size=\$(awk '/^>/{if(s)print s; s=0; next} {s+=length}END{print s}' purged.fa)
    hap_size=\$(awk '/^>/{if(s)print s; s=0; next} {s+=length}END{print s+0}' hap.fa)
    
    orig_n=\$(grep -c '^>' assembly.fa)
    purged_n=\$(grep -c '^>' purged.fa)
    hap_n=\$(grep -c '^>' hap.fa || echo 0)
    
    echo "Original assembly size: \${orig_size}" | tee -a ${haplotype_id}.purge_dups.log
    echo "Purged assembly size: \${purged_size}" | tee -a ${haplotype_id}.purge_dups.log
    echo "Haplotigs removed: \${hap_size}" | tee -a ${haplotype_id}.purge_dups.log
    echo "Original contigs: \${orig_n}" | tee -a ${haplotype_id}.purge_dups.log
    echo "Purged contigs: \${purged_n}" | tee -a ${haplotype_id}.purge_dups.log
    echo "Haplotig contigs: \${hap_n}" | tee -a ${haplotype_id}.purge_dups.log
    
    # Rename outputs to include haplotype_id
    mv purged.fa ${haplotype_id}.purged.fa
    mv hap.fa ${haplotype_id}.haplotigs.fa
    
    # Cleanup large intermediate files
    rm -f assembly.fa aligned.paf.gz split.fa split.self.paf.gz PB.stat PB.base.cov cutoffs dups.bed
    
    echo "" | tee -a ${haplotype_id}.purge_dups.log
    echo "=== purge_dups completed successfully ===" | tee -a ${haplotype_id}.purge_dups.log
    """
    
    stub:
    """
    touch ${haplotype_id}.purged.fa
    touch ${haplotype_id}.haplotigs.fa
    touch ${haplotype_id}.purge_dups.log
    """
}