/*
========================================================================================
    SNAIL PLOT MODULE (BlobToolKit assembly-stats)
========================================================================================
    Generates snail plots for assembly visualization using blobtk
    
    Input:
    - Assembly FASTA file
    - BUSCO results directory (containing full_table.tsv)
    
    Output:
    - Snail plot PNG
    - Assembly statistics JSON
    
    Can be applied to any assembly stage (contigs, scaffolds, gap-filled, etc.)
    
    Reference: https://assembly-stats.readme.io/docs/getting-started
========================================================================================
*/

process SNAIL_PLOT {
    tag "${haplotype_id}"
    label 'snail_plot'
    
    publishDir "${params.outdir}/qc/snail_plots/${qc_label}", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly_fasta), path(busco_dir), val(qc_label)
    
    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.snail.png"), emit: snail_plot
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.stats.json"), emit: stats_json
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.snail.svg"), emit: snail_svg, optional: true
    
    script:
    def busco_lineage = params.busco_lineage ?: 'actinopterygii_odb10'
    """
    set -euo pipefail
    
    # Find the BUSCO full_table.tsv file
    # BUSCO output structure: busco_dir/run_<lineage>/full_table.tsv
    BUSCO_TABLE=""
    
    # Try common BUSCO output locations
    if [ -f "${busco_dir}/run_${busco_lineage}/full_table.tsv" ]; then
        BUSCO_TABLE="${busco_dir}/run_${busco_lineage}/full_table.tsv"
    elif [ -f "${busco_dir}/full_table.tsv" ]; then
        BUSCO_TABLE="${busco_dir}/full_table.tsv"
    else
        # Search for any full_table.tsv
        BUSCO_TABLE=\$(find ${busco_dir} -name "full_table.tsv" -type f | head -n 1)
    fi
    
    if [ -z "\${BUSCO_TABLE}" ] || [ ! -f "\${BUSCO_TABLE}" ]; then
        echo "[ERROR] Could not find BUSCO full_table.tsv in ${busco_dir}" >&2
        echo "[ERROR] Directory contents:" >&2
        find ${busco_dir} -type f | head -20 >&2
        exit 1
    fi
    
    echo "[SNAIL_PLOT] Found BUSCO table: \${BUSCO_TABLE}"
    echo "[SNAIL_PLOT] Assembly: ${assembly_fasta}"
    echo "[SNAIL_PLOT] Haplotype: ${haplotype_id}"
    echo "[SNAIL_PLOT] Stage: ${qc_label}"
    
    # Run blobtk pipeline to generate snail plot
    blobtk pipeline \\
        --fasta ${assembly_fasta} \\
        --busco \${BUSCO_TABLE} \\
        --out ${haplotype_id}_${qc_label}
    
    # Rename outputs to include qc_label
    # blobtk creates: <prefix>.snail.png, <prefix>.stats.json
    if [ -f "${haplotype_id}_${qc_label}.snail.png" ]; then
        echo "[SNAIL_PLOT] Snail plot generated successfully"
    else
        echo "[ERROR] Snail plot was not generated" >&2
        ls -la >&2
        exit 1
    fi
    """
    
    stub:
    """
    touch ${haplotype_id}_${qc_label}.snail.png
    touch ${haplotype_id}_${qc_label}.stats.json
    """
}

/*
========================================================================================
    SNAIL PLOT MODULE (Minimal - Assembly Only)
========================================================================================
    Generates snail plots without BUSCO data
    Useful for quick assembly visualization
========================================================================================
*/

process SNAIL_PLOT_MINIMAL {
    tag "${haplotype_id}"
    label 'snail_plot'
    
    publishDir "${params.outdir}/qc/snail_plots/${qc_label}", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly_fasta), val(qc_label)
    
    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.snail.png"), emit: snail_plot
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.stats.json"), emit: stats_json
    
    script:
    """
    set -euo pipefail
    
    echo "[SNAIL_PLOT] Assembly: ${assembly_fasta}"
    echo "[SNAIL_PLOT] Haplotype: ${haplotype_id}"
    echo "[SNAIL_PLOT] Stage: ${qc_label}"
    
    # Run blobtk pipeline without BUSCO
    blobtk pipeline \\
        --fasta ${assembly_fasta} \\
        --out ${haplotype_id}_${qc_label}
    """
    
    stub:
    """
    touch ${haplotype_id}_${qc_label}.snail.png
    touch ${haplotype_id}_${qc_label}.stats.json
    """
}