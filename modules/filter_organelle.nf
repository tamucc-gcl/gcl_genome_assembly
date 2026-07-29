/*
========================================================================================
    FILTER ORGANELLE CONTIGS MODULE  (was FILTER_MITO_CONTIGS)
========================================================================================
    Removes organelle contigs (mitochondrial AND, for plants, plastid) from the nuclear
    assembly, using the assembled organelle sequences as bait.
    Repo location: modules/filter_organelle.nf

    Generalisation of the old mito-only filter: instead of a single mitogenome, it takes ALL
    organelle assemblies for the sample (from ORGANELLE.out.assemblies, grouped per sample),
    concatenates them into one bait, and removes any contig with high identity+coverage to it.
    This keeps organelle contigs out of purge_dups (extreme coverage), scaffolding (no Hi-C),
    and decontamination (false positives).

    Runs per-haplotype over ch_contigs (one instance handles both haplotypes).

    Input:
    - meta + assembly_fasta + organelle_fastas (one or more; concatenated into the bait)
    Output:
    - filtered:          nuclear assembly with organelle contigs removed
    - organelle_contigs: the removed contigs (for verification)
    - filter_stats:      filtering statistics
========================================================================================
*/

process FILTER_ORGANELLE {
    tag "${meta.id}"
    label 'filter_mito_contigs'

    publishDir "${params.outdir}/assembly/contig/organelle_filtered", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly_fasta), path(organelle_fastas)

    output:
    tuple val(meta), path("${meta.id}.organelle_filtered.fasta"),   emit: filtered
    tuple val(meta), path("${meta.id}.organelle_contigs.fasta"),    emit: organelle_contigs
    tuple val(meta), path("${meta.id}.organelle_filter_stats.tsv"), emit: filter_stats

    script:
    def min_identity = params.mitohifi_filter_min_identity ?: 90
    def min_coverage = params.mitohifi_filter_min_coverage ?: 50
    """
    set -euo pipefail

    echo "[FILTER_ORGANELLE] Filtering organelle contigs from ${meta.id}"
    echo "[FILTER_ORGANELLE] Min identity: ${min_identity}%  min query coverage: ${min_coverage}%"

    # Concatenate all organelle assemblies for this sample into one bait
    cat ${organelle_fastas} > organelle_bait.fasta
    echo "[FILTER_ORGANELLE] Bait sequences: \$(grep -c '^>' organelle_bait.fasta || echo 0)"

    TOTAL_CONTIGS=\$(grep -c '^>' ${assembly_fasta})
    echo "[FILTER_ORGANELLE] Total contigs in assembly: \${TOTAL_CONTIGS}"

    # 1) Align assembly contigs to the organelle bait
    minimap2 -x asm5 -t ${task.cpus} --secondary=no \\
        organelle_bait.fasta ${assembly_fasta} -o alignments.paf

    # 2) Flag contigs by identity + query coverage
    #    PAF: qname qlen qstart qend strand tname tlen tstart tend matches alnlen mapq
    awk -v min_id=${min_identity} -v min_cov=${min_coverage} '
    {
        qname=\$1; qlen=\$2; matches=\$10; alnlen=\$11
        pid  = (alnlen > 0) ? (matches/alnlen)*100 : 0
        qcov = (qlen  > 0) ? ((\$4-\$3)/qlen)*100 : 0
        if (pid >= min_id && qcov >= min_cov) {
            org[qname]=1; org_len[qname]=qlen; org_pid[qname]=pid; org_qcov[qname]=qcov
        }
    }
    END { for (c in org) print c, org_len[c], org_pid[c], org_qcov[c] }
    ' alignments.paf > organelle_contig_ids.tsv

    ORG_COUNT=\$(wc -l < organelle_contig_ids.tsv)
    echo "[FILTER_ORGANELLE] Identified \${ORG_COUNT} organelle contig(s)"

    # 3) Split into organelle and non-organelle
    if [ "\${ORG_COUNT}" -gt 0 ]; then
        cut -f1 organelle_contig_ids.tsv > organelle_names.txt
        samtools faidx ${assembly_fasta}
        xargs -a organelle_names.txt samtools faidx ${assembly_fasta} \\
            > ${meta.id}.organelle_contigs.fasta
        awk 'NR==FNR {excl[\$1]=1; next}
             /^>/ { name=substr(\$1,2); skip=(name in excl)?1:0 }
             !skip {print}
        ' organelle_names.txt ${assembly_fasta} \\
            > ${meta.id}.organelle_filtered.fasta
    else
        echo "[FILTER_ORGANELLE] No organelle contigs detected; assembly unchanged"
        cp ${assembly_fasta} ${meta.id}.organelle_filtered.fasta
        touch ${meta.id}.organelle_contigs.fasta
    fi

    # 4) Stats
    FILTERED_CONTIGS=\$(grep -c '^>' ${meta.id}.organelle_filtered.fasta || echo 0)
    ORG_TOTAL_BP=0
    [ "\${ORG_COUNT}" -gt 0 ] && ORG_TOTAL_BP=\$(grep -v '^>' ${meta.id}.organelle_contigs.fasta | tr -d '\\n' | wc -c)

    cat > ${meta.id}.organelle_filter_stats.tsv <<EOF
haplotype_id\ttotal_contigs\torganelle_contigs_removed\torganelle_total_bp\tremaining_contigs
${meta.id}\t\${TOTAL_CONTIGS}\t\${ORG_COUNT}\t\${ORG_TOTAL_BP}\t\${FILTERED_CONTIGS}
EOF
    if [ "\${ORG_COUNT}" -gt 0 ]; then
        echo "" >> ${meta.id}.organelle_filter_stats.tsv
        echo "# Removed contigs (name, length, percent_identity, query_coverage):" >> ${meta.id}.organelle_filter_stats.tsv
        cat organelle_contig_ids.tsv >> ${meta.id}.organelle_filter_stats.tsv
    fi

    echo "[FILTER_ORGANELLE] Removed \${ORG_COUNT} contig(s) (\${ORG_TOTAL_BP} bp); \${FILTERED_CONTIGS} remain"
    """

    stub:
    """
    cp ${assembly_fasta} ${meta.id}.organelle_filtered.fasta
    touch ${meta.id}.organelle_contigs.fasta
    printf 'haplotype_id\\ttotal_contigs\\torganelle_contigs_removed\\torganelle_total_bp\\tremaining_contigs\\n${meta.id}\\t100\\t0\\t0\\t100\\n' > ${meta.id}.organelle_filter_stats.tsv
    """
}
