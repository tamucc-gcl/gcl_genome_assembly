/*
========================================================================================
    GAP FILLING MODULE (TGSGapCloser)
========================================================================================
    Repo location: modules/gap_filling.nf

    Fills gaps in scaffolded assemblies using HiFi reads (TGSGapCloser, --tgstype pb, --ne).
    Input : tuple(meta, scaffold_fasta, hifi_fastq)
    Output: gap-filled assembly + detail + report (per-haplotype, meta.id).
========================================================================================
*/

process GAP_FILLING {
    tag "${meta.id}"
    label 'gap_filling'

    publishDir "${params.outdir}/assembly/scaffold/gap_filling", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(scaffold_fasta), path(hifi_fastq)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: filled_assembly
    tuple val(meta), path("${meta.id}.gap_fill_detail"), emit: details
    tuple val(meta), path("${meta.id}_gap_filling_report.txt"), emit: report

    script:
    """
    set -euo pipefail

    echo "[GAP_FILLING] Starting gap filling for ${meta.id}"
    echo "[GAP_FILLING] Started: \$(date)"
    echo "[GAP_FILLING] Input scaffold: ${scaffold_fasta}"
    echo "[GAP_FILLING] HiFi reads: ${hifi_fastq}"
    echo ""

    # Run TGSGapCloser
    tgsgapcloser \\
        --output ${meta.id} \\
        --scaff ${scaffold_fasta} \\
        --reads ${hifi_fastq} \\
        --tgstype pb \\
        --minmap_arg '-x map-hifi' \\
        --thread ${task.cpus} \\
        --ne

    echo ""
    echo "[GAP_FILLING] Gap filling complete: \$(date)"

    # Rename output file from .scaff_seqs to .fasta
    if [ -f "${meta.id}.scaff_seqs" ]; then
        mv "${meta.id}.scaff_seqs" "${meta.id}.fasta"
    else
        echo "[ERROR] Expected output file ${meta.id}.scaff_seqs not found!" >&2
        exit 1
    fi

    # Generate comprehensive report
    cat > ${meta.id}_gap_filling_report.txt <<'EOF'
================================================================================
Gap Filling Report
================================================================================
Haplotype: ${meta.id}
Generated: \$(date)

Input scaffold: ${scaffold_fasta}
HiFi reads: ${hifi_fastq}

================================================================================
INPUT ASSEMBLY STATISTICS
================================================================================
\$(if command -v seqkit &> /dev/null; then
    seqkit stats ${scaffold_fasta}
else
    echo "Number of sequences: \$(grep -c '^>' ${scaffold_fasta})"
    echo "Total length: \$(grep -v '^>' ${scaffold_fasta} | tr -d '\\n' | wc -c)"
fi)

Number of N bases (gaps) in input:
\$(grep -v '^>' ${scaffold_fasta} | tr -d '\\n' | grep -o 'N' | wc -l)

================================================================================
GAP FILLING DETAILS
================================================================================
\$(cat ${meta.id}.gap_fill_detail)

================================================================================
GAP-FILLED ASSEMBLY STATISTICS
================================================================================
\$(if command -v seqkit &> /dev/null; then
    seqkit stats ${meta.id}.fasta
else
    echo "Number of sequences: \$(grep -c '^>' ${meta.id}.fasta)"
    echo "Total length: \$(grep -v '^>' ${meta.id}.fasta | tr -d '\\n' | wc -c)"
fi)

Number of N bases (gaps) remaining:
\$(grep -v '^>' ${meta.id}.fasta | tr -d '\\n' | grep -o 'N' | wc -l)

================================================================================
GAP FILLING SUMMARY
================================================================================
Gaps before: \$(grep -v '^>' ${scaffold_fasta} | tr -d '\\n' | grep -o 'N' | wc -l)
Gaps after:  \$(grep -v '^>' ${meta.id}.fasta | tr -d '\\n' | grep -o 'N' | wc -l)
Gaps filled: \$(( \$(grep -v '^>' ${scaffold_fasta} | tr -d '\\n' | grep -o 'N' | wc -l) - \$(grep -v '^>' ${meta.id}.fasta | tr -d '\\n' | grep -o 'N' | wc -l) ))

================================================================================
OUTPUT FILES
================================================================================
Gap-filled assembly: ${meta.id}.fasta
Gap filling details: ${meta.id}.gap_fill_detail
Report: ${meta.id}_gap_filling_report.txt

================================================================================
EOF

    """

    stub:
    """
    touch ${meta.id}.fasta
    touch ${meta.id}.gap_fill_detail
    touch ${meta.id}_gap_filling_report.txt
    """
}
