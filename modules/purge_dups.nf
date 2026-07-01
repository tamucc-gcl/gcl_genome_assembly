/*
========================================================================================
    PURGE_DUPS MODULE
========================================================================================
    Remove haplotig duplications from genome assemblies using purge_dups
    Repo location: modules/purge_dups.nf

    Steps: map HiFi reads (minimap2) -> coverage (pbcstat) -> cutoffs (calcuts) ->
    self-alignment (minimap2) -> identify dups (purge_dups) -> purged assembly (get_seqs)

    Inputs:
    - meta + assembly + hifi_reads
    Outputs:
    - purged_assembly: Assembly with haplotigs removed
    - haplotigs: Removed haplotig sequences
    - log: Processing log
========================================================================================
*/

process PURGE_DUPS {
    tag "${meta.id}"
    label 'purge_dups'

    publishDir "${params.outdir}/assembly/contig/purge_dups", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly), path(hifi_reads)

    output:
    tuple val(meta), path("${meta.id}.purged.fa"), emit: purged_assembly
    tuple val(meta), path("${meta.id}.haplotigs.fa"), emit: haplotigs
    tuple val(meta), path("${meta.id}.purge_dups.log"), emit: log

    script:
    def minimap_threads = Math.max(1, task.cpus - 4)
    def pigz_threads = Math.min(4, task.cpus)
    """
    set -euo pipefail

    echo "=== purge_dups started ===" | tee ${meta.id}.purge_dups.log
    echo "Working directory: \$(pwd)" | tee -a ${meta.id}.purge_dups.log
    echo "" | tee -a ${meta.id}.purge_dups.log

    # Copy assembly to local file to avoid potential symlink issues with get_seqs
    echo "Copying assembly to local file..." | tee -a ${meta.id}.purge_dups.log
    cp ${assembly} assembly.fa

    echo "=== Step 1: Mapping HiFi reads to assembly ===" | tee -a ${meta.id}.purge_dups.log
    minimap2 -xmap-hifi -t ${minimap_threads} assembly.fa ${hifi_reads} | \\
        pigz -p ${pigz_threads} -c > aligned.paf.gz
    echo "Mapping complete" | tee -a ${meta.id}.purge_dups.log

    echo "=== Step 2: Calculating coverage statistics ===" | tee -a ${meta.id}.purge_dups.log
    pbcstat aligned.paf.gz 2>&1 | tee -a ${meta.id}.purge_dups.log

    echo "=== Step 3: Calculating coverage cutoffs ===" | tee -a ${meta.id}.purge_dups.log
    calcuts PB.stat > cutoffs 2>> ${meta.id}.purge_dups.log
    echo "Cutoffs:" | tee -a ${meta.id}.purge_dups.log
    cat cutoffs | tee -a ${meta.id}.purge_dups.log
    echo "" | tee -a ${meta.id}.purge_dups.log

    echo "=== Step 4: Splitting assembly for self-alignment ===" | tee -a ${meta.id}.purge_dups.log
    split_fa assembly.fa > split.fa
    echo "Split into \$(grep -c '^>' split.fa) sequences" | tee -a ${meta.id}.purge_dups.log

    echo "=== Step 5: Self-alignment ===" | tee -a ${meta.id}.purge_dups.log
    minimap2 -xasm5 -DP -t ${minimap_threads} split.fa split.fa | \\
        pigz -p ${pigz_threads} -c > split.self.paf.gz
    echo "Self-alignment complete" | tee -a ${meta.id}.purge_dups.log

    echo "=== Step 6: Identifying duplications ===" | tee -a ${meta.id}.purge_dups.log
    # CRITICAL: stdout goes to dups.bed, stderr (log messages) goes to log file
    # Do NOT use 2>&1 here or it corrupts the BED file!
    purge_dups -2 -T cutoffs -c PB.base.cov split.self.paf.gz > dups.bed 2>> ${meta.id}.purge_dups.log

    echo "dups.bed: \$(wc -l < dups.bed) lines" | tee -a ${meta.id}.purge_dups.log
    echo "First 10 lines of dups.bed:" | tee -a ${meta.id}.purge_dups.log
    head -10 dups.bed | tee -a ${meta.id}.purge_dups.log
    echo "" | tee -a ${meta.id}.purge_dups.log

    echo "=== Step 7: Generating purged assembly ===" | tee -a ${meta.id}.purge_dups.log

    if [ -s dups.bed ] && [ \$(wc -l < dups.bed) -gt 0 ]; then
        echo "Running: get_seqs -e dups.bed assembly.fa" | tee -a ${meta.id}.purge_dups.log
        get_seqs -e dups.bed assembly.fa 2>&1 | tee -a ${meta.id}.purge_dups.log
        echo "get_seqs completed" | tee -a ${meta.id}.purge_dups.log
    else
        echo "No duplications identified, using original assembly" | tee -a ${meta.id}.purge_dups.log
        cp assembly.fa purged.fa
        touch hap.fa
    fi

    # Calculate summary statistics
    echo "" | tee -a ${meta.id}.purge_dups.log
    echo "=== Summary ===" | tee -a ${meta.id}.purge_dups.log

    orig_size=\$(awk '/^>/{if(s)print s; s=0; next} {s+=length}END{print s}' assembly.fa)
    purged_size=\$(awk '/^>/{if(s)print s; s=0; next} {s+=length}END{print s}' purged.fa)
    hap_size=\$(awk '/^>/{if(s)print s; s=0; next} {s+=length}END{print s+0}' hap.fa)

    orig_n=\$(grep -c '^>' assembly.fa)
    purged_n=\$(grep -c '^>' purged.fa)
    hap_n=\$(grep -c '^>' hap.fa || echo 0)

    echo "Original assembly size: \${orig_size}" | tee -a ${meta.id}.purge_dups.log
    echo "Purged assembly size: \${purged_size}" | tee -a ${meta.id}.purge_dups.log
    echo "Haplotigs removed: \${hap_size}" | tee -a ${meta.id}.purge_dups.log
    echo "Original contigs: \${orig_n}" | tee -a ${meta.id}.purge_dups.log
    echo "Purged contigs: \${purged_n}" | tee -a ${meta.id}.purge_dups.log
    echo "Haplotig contigs: \${hap_n}" | tee -a ${meta.id}.purge_dups.log

    # Rename outputs to include haplotype id
    mv purged.fa ${meta.id}.purged.fa
    mv hap.fa ${meta.id}.haplotigs.fa

    # Cleanup large intermediate files
    rm -f assembly.fa aligned.paf.gz split.fa split.self.paf.gz PB.stat PB.base.cov cutoffs dups.bed

    echo "" | tee -a ${meta.id}.purge_dups.log
    echo "=== purge_dups completed successfully ===" | tee -a ${meta.id}.purge_dups.log
    """

    stub:
    """
    touch ${meta.id}.purged.fa
    touch ${meta.id}.haplotigs.fa
    touch ${meta.id}.purge_dups.log
    """
}
