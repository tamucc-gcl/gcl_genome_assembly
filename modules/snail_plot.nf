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
    
    publishDir "${params.outdir}/qc/snail_plots/${qc_label}", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly_fasta), path(busco_dir), val(qc_label)
    
    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.fasta"), emit: assembly
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_busco"), emit: busco
    
    script:
    """
    # Copy assembly with standardized name
    cp ${assembly_fasta} ${haplotype_id}_${qc_label}.fasta
    
    # Copy BUSCO directory with standardized name
    cp -r ${busco_dir} ${haplotype_id}_${qc_label}_busco
    
    echo "[SNAIL_PLOT] Staged files for ${haplotype_id} (${qc_label})"
    echo "[SNAIL_PLOT] Assembly: ${haplotype_id}_${qc_label}.fasta"
    echo "[SNAIL_PLOT] BUSCO dir: ${haplotype_id}_${qc_label}_busco"
    """
    
    stub:
    """
    touch ${haplotype_id}_${qc_label}.fasta
    mkdir -p ${haplotype_id}_${qc_label}_busco
    """
}

/*
========================================================================================
    SNAIL PLOT MODULE (Minimal - Assembly Only)
========================================================================================
    Stages assembly for snail plot generation without BUSCO data
========================================================================================
*/

process SNAIL_PLOT_MINIMAL {
    tag "${haplotype_id}"
    label 'snail_plot'
    
    publishDir "${params.outdir}/qc/snail_plots/${qc_label}", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly_fasta), val(qc_label)
    
    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.fasta"), emit: assembly
    
    script:
    """
    cp ${assembly_fasta} ${haplotype_id}_${qc_label}.fasta
    
    echo "[SNAIL_PLOT] Staged assembly for ${haplotype_id} (${qc_label})"
    """
    
    stub:
    """
    touch ${haplotype_id}_${qc_label}.fasta
    """
}