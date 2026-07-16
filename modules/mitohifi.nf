/*
========================================================================================
    MITOHIFI MODULE
========================================================================================
    Assembles and annotates the mitochondrial genome from PacBio HiFi reads
    using MitoHiFi in reads mode (-r).
    Repo location: modules/mitohifi.nf

    Runs concurrently with HIFIASM — both consume HiFi reads independently.
    The assembled mitogenome is used downstream by FILTER_MITO_CONTIGS to
    remove mitochondrial sequences from the nuclear assembly.

    Input:
    - meta + hifi_fastq: HiFi reads (from BAM_TO_FASTQ)
    - ref_fasta: Closely related mitogenome FASTA (from FIND_MITO_REFERENCE)
    - ref_gb: Matching GenBank annotation (from FIND_MITO_REFERENCE)

    Output:
    - mitogenome: Final assembled + rotated mitogenome FASTA
    - annotation: GenBank annotation of the mitogenome
    - stats: Assembly statistics and QC metrics
    - contigs_fasta: All candidate mitochondrial contigs (for filtering)
========================================================================================
*/

process MITOHIFI {
    tag "${meta.sample}"
    label 'mitohifi'

    publishDir "${params.outdir}/mitogenome", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hifi_fastq), path(ref_fasta), path(ref_gb)

    output:
    tuple val(meta), path("${meta.sample}_mitogenome.fasta"),           emit: mitogenome
    tuple val(meta), path("${meta.sample}_mitogenome.gb"),              emit: annotation
    tuple val(meta), path("${meta.sample}_mito_stats.tsv"),             emit: stats
    tuple val(meta), path("${meta.sample}_mito_contigs.fasta"),         emit: contigs_fasta
    //tuple val(meta), path("mitohifi_output"),                          emit: output_dir
    tuple val(meta), path("${meta.sample}_final_mitogenome*.png"),      emit: gene_map, optional: true

    script:
    def genetic_code = params.mitohifi_genetic_code ?: 2
    def perc_id      = params.mitohifi_perc_identity ?: 50
    def cov_cutoff   = params.mitohifi_cov_cutoff ?: 'auto'
    def bloom_filter = params.mitohifi_bloom_filter ? '--bloom-filter' : ''

    // Build optional arguments
    def cov_arg = cov_cutoff == 'auto' ? '' : "--mitos-cov ${cov_cutoff}"

    """
    set -euo pipefail

    echo "[MITOHIFI] Assembling mitogenome for ${meta.sample}"
    echo "[MITOHIFI] Genetic code: ${genetic_code}"
    echo "[MITOHIFI] Min percent identity: ${perc_id}%"
    echo "[MITOHIFI] Reference: ${ref_fasta}"
    echo "[MITOHIFI] Started: \$(date)"

    # Create output directory
    mkdir -p mitohifi_output

    # Run MitoHiFi in reads mode
    mitohifi.py \\
        -r ${hifi_fastq} \\
        -f ${ref_fasta} \\
        -g ${ref_gb} \\
        -o ${genetic_code} \\
        -t ${task.cpus} \\
        -p ${perc_id} \\
        ${bloom_filter} \\
        ${cov_arg}

    echo "[MITOHIFI] MitoHiFi complete: \$(date)"

    # -------------------------------------------------------------------------
    # Collect and rename outputs with sample prefix
    # -------------------------------------------------------------------------

    # Final mitogenome
    if [ -f final_mitogenome.fasta ]; then
        # Prefix the header with sample for clarity
        sed "s/^>/>mitogenome_${meta.sample} /" final_mitogenome.fasta \\
            > ${meta.sample}_mitogenome.fasta
    else
        echo "[MITOHIFI] WARNING: final_mitogenome.fasta not found"
        echo "[MITOHIFI] MitoHiFi may have failed to assemble a complete mitogenome"
        echo "[MITOHIFI] Checking for partial results..."
        ls -la *.fasta 2>/dev/null || true
        exit 1
    fi

    # GenBank annotation
    if [ -f final_mitogenome.gb ]; then
        cp final_mitogenome.gb ${meta.sample}_mitogenome.gb
    else
        touch ${meta.sample}_mitogenome.gb
    fi

    # Gene map image(s)
    for png in final_mitogenome*.png; do
        [ -f "\${png}" ] && cp "\${png}" "${meta.sample}_\${png}"
    done

    # All candidate mito contigs (used by FILTER_MITO_CONTIGS)
    if [ -f contigs_filtering/all_potential_contigs.fa ]; then
        cp contigs_filtering/all_potential_contigs.fa ${meta.sample}_mito_contigs.fasta
    elif [ -f final_mitogenome.fasta ]; then
        # Fallback: use the final mitogenome itself as the filter target
        cp final_mitogenome.fasta ${meta.sample}_mito_contigs.fasta
    fi

    # -------------------------------------------------------------------------
    # Compile statistics
    # -------------------------------------------------------------------------
    MITO_LEN=\$(grep -v '^>' ${meta.sample}_mitogenome.fasta | tr -d '\\n' | wc -c)

    # Check circularity from MitoHiFi log
    CIRCULAR="unknown"
    if grep -qi "circular" contigs_stats.tsv 2>/dev/null; then
        CIRCULAR="yes"
    fi

    # Count annotated genes
    GENE_COUNT=0
    if [ -f ${meta.sample}_mitogenome.gb ]; then
        GENE_COUNT=\$(grep -c '^ */gene=' ${meta.sample}_mitogenome.gb 2>/dev/null || echo 0)
    fi

    # Count tRNAs
    TRNA_COUNT=0
    if [ -f ${meta.sample}_mitogenome.gb ]; then
        TRNA_COUNT=\$(grep -c '^ */product="tRNA' ${meta.sample}_mitogenome.gb 2>/dev/null || echo 0)
    fi

    # Count rRNAs
    RRNA_COUNT=0
    if [ -f ${meta.sample}_mitogenome.gb ]; then
        RRNA_COUNT=\$(grep -c 'rRNA' ${meta.sample}_mitogenome.gb 2>/dev/null || echo 0)
    fi

    cat > ${meta.sample}_mito_stats.tsv <<EOF
sample_id\tmitogenome_length\tcircular\tgene_count\ttrna_count\trrna_count\tgenetic_code
${meta.sample}\t\${MITO_LEN}\t\${CIRCULAR}\t\${GENE_COUNT}\t\${TRNA_COUNT}\t\${RRNA_COUNT}\t${genetic_code}
EOF

    echo "[MITOHIFI] Mitogenome summary:"
    echo "  Length:   \${MITO_LEN} bp"
    echo "  Circular: \${CIRCULAR}"
    echo "  Genes:    \${GENE_COUNT} (tRNAs: \${TRNA_COUNT}, rRNAs: \${RRNA_COUNT})"

    # Move all MitoHiFi intermediates into output directory
    for f in contigs_stats.tsv shared_genes.tsv contigs_filtering/ reads_mapping_and_assembly/; do
        [ -e "\${f}" ] && mv "\${f}" mitohifi_output/ 2>/dev/null || true
    done

    echo "[MITOHIFI] Complete: \$(date)"
    """

    stub:
    """
    mkdir -p mitohifi_output
    echo ">mitogenome_${meta.sample} stub" > ${meta.sample}_mitogenome.fasta
    echo "ATCGATCG" >> ${meta.sample}_mitogenome.fasta
    touch ${meta.sample}_mitogenome.gb
    touch ${meta.sample}_mito_stats.tsv
    cp ${meta.sample}_mitogenome.fasta ${meta.sample}_mito_contigs.fasta
    touch ${meta.sample}_final_mitogenome.annotation.png
    """
}
