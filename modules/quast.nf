/*
========================================================================================
    QUAST MODULE
========================================================================================
    Assembly quality assessment for both haplotypes
========================================================================================
*/

process QUAST {
    tag "${sample_id}"
    label 'quast'
    
    publishDir "${params.outdir}/${sample_id}/qc/assembly/quast", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hap1_fasta), path(hap2_fasta)
    
    output:
    tuple val(sample_id), path("quast_results"), emit: results
    
    script:
    """
    quast.py \\
        ${hap1_fasta} \\
        ${hap2_fasta} \\
        --threads ${task.cpus} \\
        --labels ${sample_id}.hap1,${sample_id}.hap2 \\
        --output-dir quast_results
    """
    
    stub:
    """
    mkdir -p quast_results
    touch quast_results/report.txt
    touch quast_results/report.html
    """
}