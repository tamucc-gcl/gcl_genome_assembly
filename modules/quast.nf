/*
========================================================================================
    QUAST MODULE
========================================================================================
    Assembly quality assessment for both haplotypes
========================================================================================
*/

process QUAST {
    tag "${sample_id}_${qc_label}"
    label 'quast'
    
    //publishDir "${params.outdir}/qc/assembly/${qc_label}/quast", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hap1_fasta), path(hap2_fasta), val(qc_label)
    
    output:
    tuple val(sample_id), path("${sample_id}_quast"), val(qc_label), emit: results
    
    script:
    """
    quast.py \\
        ${hap1_fasta} \\
        ${hap2_fasta} \\
        --threads ${task.cpus} \\
        --labels ${sample_id}.hap1,${sample_id}.hap2 \\
        --output-dir ${sample_id}_quast
    """
    
    stub:
    """
    mkdir -p ${sample_id}_quast
    touch ${sample_id}_quast/report.txt
    touch ${sample_id}_quast/report.html
    """
}