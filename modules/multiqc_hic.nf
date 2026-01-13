/*
========================================================================================
    MULTIQC HI-C MODULE
========================================================================================
    Aggregates FastQC reports from Hi-C reads into a single MultiQC report
========================================================================================
*/

process MULTIQC_HIC {
    label 'multiqc'
    
    publishDir "${params.outdir}/qc/hic/${qc_label}", mode: params.publish_dir_mode
    
    input:
    path(fastqc_zips)  // all FastQC zip files collected
    val qc_label       // "raw" or "trimmed"
    
    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"),        emit: data
    
    script:
    def title = qc_label == "raw" ? "Hi-C Raw Reads QC" : "Hi-C Trimmed Reads QC"
    def comment = qc_label == "raw" ? 
        "Quality control metrics for Hi-C paired-end reads before trimming" : 
        "Quality control metrics for Hi-C paired-end reads after trimming"
    """
    multiqc \\
        --force \\
        --filename multiqc_report.html \\
        --title "${title}" \\
        --comment "${comment}" \\
        .
    """
    
    stub:
    """
    mkdir -p multiqc_report_data
    touch multiqc_report.html
    touch multiqc_report_data/multiqc_data.json
    """
}