/*
========================================================================================
    FILTER MITO CONTIGS MODULE
========================================================================================
    Removes mitochondrial contigs from the nuclear assembly using the MitoHiFi
    assembled mitogenome as reference.
    Repo location: modules/filter_mito_contigs.nf

    Uses minimap2 to align contigs against the mitogenome, then filters out any
    contig with high identity and coverage alignment to the mitogenome. This
    prevents mito contigs from confusing purge_dups (extreme coverage),
    scaffolding (no Hi-C signal), and decontamination (false positive flags).

    Runs per-haplotype over ch_contigs (one instance handles both haplotypes;
    no HAP1/HAP2 aliasing needed now that the assembly is forked upstream).

    Input:
    - meta + assembly_fasta + mitogenome_fasta
    Output:
    - filtered: Nuclear assembly with mito contigs removed
    - mito_contigs: Extracted mitochondrial contigs (for verification)
    - filter_stats: Filtering statistics (contigs removed, sizes, identities)
========================================================================================
*/

process FILTER_MITO_CONTIGS {
    tag "${meta.id}"
    label 'filter_mito_contigs'

    publishDir "${params.outdir}/assembly/contig/mito_filtered", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly_fasta), path(mitogenome_fasta)

    output:
    tuple val(meta), path("${meta.id}.mito_filtered.fasta"),  emit: filtered
    tuple val(meta), path("${meta.id}.mito_contigs.fasta"),   emit: mito_contigs
    tuple val(meta), path("${meta.id}.mito_filter_stats.tsv"),emit: filter_stats

    script:
    def min_identity = params.mitohifi_filter_min_identity ?: 90
    def min_coverage = params.mitohifi_filter_min_coverage ?: 50
    """
    set -euo pipefail

    echo "[FILTER_MITO] Filtering mitochondrial contigs from ${meta.id}"
    echo "[FILTER_MITO] Min identity: ${min_identity}%"
    echo "[FILTER_MITO] Min query coverage: ${min_coverage}%"
    echo "[FILTER_MITO] Started: \$(date)"

    TOTAL_CONTIGS=\$(grep -c '^>' ${assembly_fasta})
    echo "[FILTER_MITO] Total contigs in assembly: \${TOTAL_CONTIGS}"

    # -------------------------------------------------------------------------
    # 1) Align assembly contigs to mitogenome
    # -------------------------------------------------------------------------
    minimap2 \\
        -x asm5 \\
        -t ${task.cpus} \\
        --secondary=no \\
        ${mitogenome_fasta} \\
        ${assembly_fasta} \\
        -o alignments.paf

    # -------------------------------------------------------------------------
    # 2) Identify mito contigs based on identity and coverage thresholds
    #    PAF columns: qname qlen qstart qend strand tname tlen tstart tend matches alnlen mapq
    # -------------------------------------------------------------------------
    awk -v min_id=${min_identity} -v min_cov=${min_coverage} '
    BEGIN { OFS="\\t" }
    {
        qname  = \$1
        qlen   = \$2
        qstart = \$3
        qend   = \$4
        matches = \$10
        alnlen  = \$11

        # Percent identity
        pid = (alnlen > 0) ? (matches / alnlen) * 100 : 0

        # Query coverage (how much of the contig aligns to the mitogenome)
        qcov = (qlen > 0) ? ((\$4 - \$3) / qlen) * 100 : 0

        # Accumulate per-contig (in case of multiple alignments)
        if (pid >= min_id && qcov >= min_cov) {
            mito[qname] = 1
            mito_pid[qname] = pid
            mito_qcov[qname] = qcov
            mito_qlen[qname] = qlen
        }
    }
    END {
        for (c in mito) {
            print c, mito_qlen[c], mito_pid[c], mito_qcov[c]
        }
    }
    ' alignments.paf > mito_contig_ids.tsv

    MITO_COUNT=\$(wc -l < mito_contig_ids.tsv)
    echo "[FILTER_MITO] Identified \${MITO_COUNT} mitochondrial contig(s)"

    # -------------------------------------------------------------------------
    # 3) Extract mito and non-mito contigs
    # -------------------------------------------------------------------------
    if [ "\${MITO_COUNT}" -gt 0 ]; then
        # Get contig names to exclude
        cut -f1 mito_contig_ids.tsv > mito_names.txt

        # Extract mito contigs
        samtools faidx ${assembly_fasta}
        xargs -a mito_names.txt samtools faidx ${assembly_fasta} \\
            > ${meta.id}.mito_contigs.fasta

        # Build filtered assembly (all contigs NOT in the mito list)
        awk 'NR==FNR {exclude[\$1]=1; next}
             /^>/ {
                 name = substr(\$1, 2)
                 skip = (name in exclude) ? 1 : 0
             }
             !skip {print}
        ' mito_names.txt ${assembly_fasta} \\
            > ${meta.id}.mito_filtered.fasta
    else
        # No mito contigs found — pass through unchanged
        echo "[FILTER_MITO] No mitochondrial contigs detected; assembly unchanged"
        cp ${assembly_fasta} ${meta.id}.mito_filtered.fasta
        touch ${meta.id}.mito_contigs.fasta
    fi

    # -------------------------------------------------------------------------
    # 4) Compile statistics
    # -------------------------------------------------------------------------
    FILTERED_CONTIGS=\$(grep -c '^>' ${meta.id}.mito_filtered.fasta || echo 0)
    MITO_TOTAL_BP=0
    if [ "\${MITO_COUNT}" -gt 0 ]; then
        MITO_TOTAL_BP=\$(grep -v '^>' ${meta.id}.mito_contigs.fasta | tr -d '\\n' | wc -c)
    fi

    cat > ${meta.id}.mito_filter_stats.tsv <<EOF
haplotype_id\ttotal_contigs\tmito_contigs_removed\tmito_total_bp\tremaining_contigs
${meta.id}\t\${TOTAL_CONTIGS}\t\${MITO_COUNT}\t\${MITO_TOTAL_BP}\t\${FILTERED_CONTIGS}
EOF

    # Append per-contig details if any mito contigs were found
    if [ "\${MITO_COUNT}" -gt 0 ]; then
        echo "" >> ${meta.id}.mito_filter_stats.tsv
        echo "# Removed contigs (contig_name, length, percent_identity, query_coverage):" \\
            >> ${meta.id}.mito_filter_stats.tsv
        cat mito_contig_ids.tsv >> ${meta.id}.mito_filter_stats.tsv
    fi

    echo "[FILTER_MITO] Filtering complete:"
    echo "  Input contigs:   \${TOTAL_CONTIGS}"
    echo "  Mito removed:    \${MITO_COUNT} (\${MITO_TOTAL_BP} bp)"
    echo "  Output contigs:  \${FILTERED_CONTIGS}"
    echo "[FILTER_MITO] Complete: \$(date)"
    """

    stub:
    """
    cp ${assembly_fasta} ${meta.id}.mito_filtered.fasta
    touch ${meta.id}.mito_contigs.fasta
    printf 'haplotype_id\\ttotal_contigs\\tmito_contigs_removed\\tmito_total_bp\\tremaining_contigs\\n${meta.id}\\t100\\t0\\t0\\t100\\n' > ${meta.id}.mito_filter_stats.tsv
    """
}
