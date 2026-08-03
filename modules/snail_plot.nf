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
    tuple val(haplotype_id), path(assembly_fasta), path(busco_dir), val(qc_label)
    
    output:
    //tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.fasta"), emit: assembly
    //tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_busco"), emit: busco
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_snail.svg"), emit: snail
    
    script:
    """
    set -euo pipefail

    # Create the BlobDir from the assembly
    blobtools create \\
        --fasta ${assembly_fasta} \\
        ${haplotype_id}_${qc_label}

    # BUSCO lineage is taxonomy-derived per sample, so the run_<lineage> dir name varies --
    # discover the full_table.tsv rather than assuming it. -print -quit grabs the first match
    # without a head pipe (avoids SIGPIPE under pipefail).
    FULL_TABLE=\$(find ${busco_dir} -name full_table.tsv -print -quit)
    if [ -z "\${FULL_TABLE}" ]; then
        echo "[SNAIL_PLOT] ERROR: no full_table.tsv found under ${busco_dir}" >&2
        exit 1
    fi

    blobtools add \\
        --busco "\${FULL_TABLE}" \\
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
