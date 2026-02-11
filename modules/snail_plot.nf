/*
========================================================================================
    SNAIL PLOT MODULE (BlobToolKit)
========================================================================================
    Generates snail plots for assembly visualization using blobtools2
    
    Input:
    - Assembly FASTA file
    - BUSCO results directory (containing full_table.tsv)
    
    Output:
    - Snail plot PNG
    - Assembly statistics JSON
    
    Can be applied to any assembly stage (contigs, scaffolds, gap-filled, etc.)
    
    Reference: https://blobtoolkit.genomehubs.org/blobtools2/
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
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_blobdir"), emit: blobdir
    
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
    
    # Create BlobDir with assembly
    blobtools create \\
        --fasta ${assembly_fasta} \\
        ${haplotype_id}_${qc_label}_blobdir
    
    # Add BUSCO data
    blobtools add \\
        --busco \${BUSCO_TABLE} \\
        ${haplotype_id}_${qc_label}_blobdir
    
    # Generate snail plot
    # blobtools view generates static images
    blobtools view \\
        --plot \\
        --view snail \\
        --format png \\
        --out . \\
        ${haplotype_id}_${qc_label}_blobdir
    
    # Find and rename the generated snail plot
    # blobtools view creates files like: <blobdir_name>.snail.png or snail.png
    SNAIL_FILE=""
    if [ -f "${haplotype_id}_${qc_label}_blobdir.snail.png" ]; then
        SNAIL_FILE="${haplotype_id}_${qc_label}_blobdir.snail.png"
    elif [ -f "snail.png" ]; then
        SNAIL_FILE="snail.png"
    else
        # Search for any snail.png
        SNAIL_FILE=\$(find . -maxdepth 1 -name "*snail*.png" -type f | head -n 1)
    fi
    
    if [ -z "\${SNAIL_FILE}" ] || [ ! -f "\${SNAIL_FILE}" ]; then
        echo "[ERROR] Snail plot was not generated" >&2
        echo "[ERROR] Directory contents:" >&2
        ls -la >&2
        exit 1
    fi
    
    # Rename to standard output name
    mv "\${SNAIL_FILE}" "${haplotype_id}_${qc_label}.snail.png"
    
    echo "[SNAIL_PLOT] Snail plot generated successfully"
    """
    
    stub:
    """
    mkdir -p ${haplotype_id}_${qc_label}_blobdir
    touch ${haplotype_id}_${qc_label}.snail.png
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
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_blobdir"), emit: blobdir
    
    script:
    """
    set -euo pipefail
    
    echo "[SNAIL_PLOT] Assembly: ${assembly_fasta}"
    echo "[SNAIL_PLOT] Haplotype: ${haplotype_id}"
    echo "[SNAIL_PLOT] Stage: ${qc_label}"
    
    # Create BlobDir with assembly only
    blobtools create \\
        --fasta ${assembly_fasta} \\
        ${haplotype_id}_${qc_label}_blobdir
    
    # Generate snail plot
    blobtools view \\
        --plot \\
        --view snail \\
        --format png \\
        --out . \\
        ${haplotype_id}_${qc_label}_blobdir
    
    # Find and rename the generated snail plot
    SNAIL_FILE=""
    if [ -f "${haplotype_id}_${qc_label}_blobdir.snail.png" ]; then
        SNAIL_FILE="${haplotype_id}_${qc_label}_blobdir.snail.png"
    elif [ -f "snail.png" ]; then
        SNAIL_FILE="snail.png"
    else
        SNAIL_FILE=\$(find . -maxdepth 1 -name "*snail*.png" -type f | head -n 1)
    fi
    
    if [ -z "\${SNAIL_FILE}" ] || [ ! -f "\${SNAIL_FILE}" ]; then
        echo "[ERROR] Snail plot was not generated" >&2
        ls -la >&2
        exit 1
    fi
    
    mv "\${SNAIL_FILE}" "${haplotype_id}_${qc_label}.snail.png"
    
    echo "[SNAIL_PLOT] Snail plot generated successfully"
    """
    
    stub:
    """
    mkdir -p ${haplotype_id}_${qc_label}_blobdir
    touch ${haplotype_id}_${qc_label}.snail.png
    """
}