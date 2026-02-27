/*
========================================================================================
    GAP FILLING MODULE (TGSGapCloser)
========================================================================================
    Fills gaps in scaffolded assemblies using HiFi reads
    
    Input:
    - Scaffolded assembly FASTA (from YaHS round 1 or 2)
    - HiFi reads (from BAM_TO_FASTQ)
    
    Output:
    - Gap-filled assembly: ${haplotype_id}.fasta
    - Gap filling details: ${haplotype_id}.gap_fill_detail
    
    Notes:
    - Uses TGSGapCloser with PacBio HiFi reads (--tgstype pb)
    - Runs with --ne flag (no error correction, appropriate for HiFi)
    - Works on final scaffolded assemblies from any pipeline stage
========================================================================================
*/

process GAP_FILLING {
    tag "${haplotype_id}"
    label 'gap_filling'
    
    publishDir "${params.outdir}/assembly/final", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(scaffold_fasta), path(hifi_fastq)
    
    output:
    tuple val(haplotype_id), path("${haplotype_id}.fasta"), emit: filled_assembly
    tuple val(haplotype_id), path("${haplotype_id}.gap_fill_detail"), emit: details
    tuple val(haplotype_id), path("${haplotype_id}_gap_filling_report.txt"), emit: report
    
    script:
    """
    set -euo pipefail
    
    echo "[GAP_FILLING] Starting gap filling for ${haplotype_id}"
    echo "[GAP_FILLING] Started: \$(date)"
    echo "[GAP_FILLING] Input scaffold: ${scaffold_fasta}"
    echo "[GAP_FILLING] HiFi reads: ${hifi_fastq}"
    echo ""
    
    # Run TGSGapCloser
    tgsgapcloser \\
        --output ${haplotype_id} \\
        --scaff ${scaffold_fasta} \\
        --reads ${hifi_fastq} \\
        --tgstype pb \\
        --minmap_arg '-x map-hifi' \\
        --thread ${task.cpus} \\
        --ne
    
    echo ""
    echo "[GAP_FILLING] Gap filling complete: \$(date)"
    
    # Rename output file from .scaff_seqs to .fasta
    if [ -f "${haplotype_id}.scaff_seqs" ]; then
        mv "${haplotype_id}.scaff_seqs" "${haplotype_id}.fasta"
    else
        echo "[ERROR] Expected output file ${haplotype_id}.scaff_seqs not found!" >&2
        exit 1
    fi
    
    # Generate comprehensive report
    cat > ${haplotype_id}_gap_filling_report.txt <<'EOF'
================================================================================
Gap Filling Report
================================================================================
Haplotype: ${haplotype_id}
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
\$(cat ${haplotype_id}.gap_fill_detail)

================================================================================
GAP-FILLED ASSEMBLY STATISTICS
================================================================================
\$(if command -v seqkit &> /dev/null; then
    seqkit stats ${haplotype_id}.fasta
else
    echo "Number of sequences: \$(grep -c '^>' ${haplotype_id}.fasta)"
    echo "Total length: \$(grep -v '^>' ${haplotype_id}.fasta | tr -d '\\n' | wc -c)"
fi)

Number of N bases (gaps) remaining:
\$(grep -v '^>' ${haplotype_id}.fasta | tr -d '\\n' | grep -o 'N' | wc -l)

================================================================================
GAP FILLING SUMMARY
================================================================================
Gaps before: \$(grep -v '^>' ${scaffold_fasta} | tr -d '\\n' | grep -o 'N' | wc -l)
Gaps after:  \$(grep -v '^>' ${haplotype_id}.fasta | tr -d '\\n' | grep -o 'N' | wc -l)
Gaps filled: \$(( \$(grep -v '^>' ${scaffold_fasta} | tr -d '\\n' | grep -o 'N' | wc -l) - \$(grep -v '^>' ${haplotype_id}.fasta | tr -d '\\n' | grep -o 'N' | wc -l) ))

================================================================================
OUTPUT FILES
================================================================================
Gap-filled assembly: ${haplotype_id}.fasta
Gap filling details: ${haplotype_id}.gap_fill_detail
Report: ${haplotype_id}_gap_filling_report.txt

================================================================================
EOF
    
    """
    
    stub:
    """
    touch ${haplotype_id}.fasta
    touch ${haplotype_id}.gap_fill_detail
    touch ${haplotype_id}_gap_filling_report.txt
    """
}