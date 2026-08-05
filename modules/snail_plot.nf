/*
========================================================================================
    SNAIL PLOT MODULE (BlobToolKit)
========================================================================================
    Generates snail plots for assembly visualization using blobtools2
    
    Input:
    - Assembly FASTA file
    - BUSCO results directory (containing full_table.tsv)
    
    Output:
    - Staged assembly FASTA
    - Staged BUSCO directory
    
    Can be applied to any assembly stage (contigs, scaffolds, gap-filled, etc.)
    
    TODO: Add actual snail plot generation once command is debugged
========================================================================================
*/

process SNAIL_PLOT {
    tag "${haplotype_id}"
    label 'snail_plot'
    
    publishDir "${params.outdir}/snail_plots", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly_fasta), path(busco_full_table), val(qc_label)

    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_snail.svg"), emit: snail

    script:
    """
    set -euo pipefail

    blobtools create \\
        --fasta ${assembly_fasta} \\
        ${haplotype_id}_${qc_label}

    blobtools add \\
        --busco ${busco_full_table} \\
        ${haplotype_id}_${qc_label}

    blobtk plot \\
        --blobdir ${haplotype_id}_${qc_label} \\
        --view snail \\
        --output ${haplotype_id}_${qc_label}_snail.svg
    """
    
    stub:
    """
    touch ${haplotype_id}_${qc_label}_snail.svg
    """
}
