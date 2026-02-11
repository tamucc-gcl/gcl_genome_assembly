/*
========================================================================================
    TELOMERE SCANNING MODULE
========================================================================================
    Scans scaffold ends for telomere motif repeats.
    
    Input : tuple(haplotype_id, scaffold_fasta)
    Output: per-scaffold TSV with telomere presence at 5' and 3' ends
    
    Parameters (flattened):
      params.telomere_motif       - Base telomere motif (default: TTAGGG)
                                    Complement/reverse variants generated automatically
      params.telomere_window      - Window size in bp to search at each end (default: 10000)
      params.telomere_min_repeats - Minimum consecutive repeats required (default: 10)
========================================================================================
*/

process SCAN_TELOMERES {
    tag "${haplotype_id}"
    label 'telomere_scan'

    //publishDir "${params.outdir}/qc/telomeres", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(scaffold_fasta)

    output:
    tuple val(haplotype_id), path("${haplotype_id}.telomeres.tsv"), emit: telomeres
    tuple val(haplotype_id), path("${haplotype_id}.telomere_summary.tsv"), emit: summary

    script:
    def motif = params.telomere_motif
    def window = params.telomere_window
    def min_repeats = params.telomere_min_repeats
    """
    set -euo pipefail

    ${projectDir}/py_scripts/scan_telomeres.py \\
        ${scaffold_fasta} \\
        --id ${haplotype_id} \\
        --motif ${motif} \\
        --window ${window} \\
        --min-repeats ${min_repeats} \\
        --output ${haplotype_id}.telomeres.tsv

    # Generate summary statistics
    awk -F'\\t' '
        BEGIN {
            OFS="\\t"
        }
        NR == 1 { next }
        {
            total++
            total_len += \$3
            if (\$4 == "1") five++
            if (\$5 == "1") three++
            if (\$6 == "1") both++
        }
        END {
            five_only = five - both
            three_only = three - both
            none = total - five - three + both
            print "haplotype_id", "scaffolds", "total_length", "telomere_5prime", "telomere_3prime", "telomere_both", "telomere_5prime_only", "telomere_3prime_only", "telomere_none", "pct_with_telomere"
            pct = (total > 0) ? sprintf("%.2f", (five + three - both) / total * 100) : "0.00"
            print "${haplotype_id}", total, total_len, five+0, three+0, both+0, five_only+0, three_only+0, none+0, pct
        }
    ' ${haplotype_id}.telomeres.tsv > ${haplotype_id}.telomere_summary.tsv
    """
}