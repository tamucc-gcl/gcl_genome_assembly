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
    
    //publishDir "${params.outdir}/qc/assembly/quast", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hap1_fasta), path(hap2_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_quast"), emit: results
    
    script:
    """
    quast.py \\
        ${hap1_fasta} \\
        ${hap2_fasta} \\
        --threads ${task.cpus} \\
        --labels ${sample_id}.hap1,${sample_id}.hap2 \\
        --output-dir ${sample_id}_quast \\
        --contig-thresholds 10000,50000,100000,500000,1000000,5000000,10000000
    """
    
    stub:
    """
    mkdir -p ${sample_id}_quast
    touch ${sample_id}_quast/report.txt
    touch ${sample_id}_quast/report.html
    """
}

/*
========================================================================================
    QUAST FINAL MODULE
========================================================================================
    Assembly quality assessment for all final gap-filled genomes
    - Compares all samples and haplotypes in a single report
    - Useful for cross-sample comparison and final QC summary
========================================================================================
*/

process QUAST_FINAL {
    tag "all_samples"
    label 'quast'
    
    publishDir "${params.outdir}/qc/assembly", mode: params.publish_dir_mode
    
    input:
    path(assemblies)  // All gap-filled assemblies collected
    val(labels)       // Corresponding labels (sample_id_hap1, sample_id_hap2, etc.)
    
    output:
    path("quast_final"), emit: results
    path("quast_final/report.tsv"), emit: report_tsv
    path("quast_final/report.html"), emit: report_html
    
    script:
    def labels_str = labels.join(',')
    """
    quast.py \\
        ${assemblies} \\
        --threads ${task.cpus} \\
        --labels ${labels_str} \\
        --output-dir quast_final \\
        --contig-thresholds 10000,50000,100000,500000,1000000,5000000,10000000
    """
    
    stub:
    """
    mkdir -p quast_final
    touch quast_final/report.txt
    touch quast_final/report.tsv
    touch quast_final/report.html
    """
}