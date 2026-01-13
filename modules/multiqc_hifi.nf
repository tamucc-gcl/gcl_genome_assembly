/*
========================================================================================
    MULTIQC HiFi MODULE
========================================================================================
    Aggregates FastQC reports from Hi-C reads into a single MultiQC report
========================================================================================
*/

process MULTIQC_HIFI {
    label 'multiqc'
    
    publishDir "${params.outdir}/qc/hifi", mode: params.publish_dir_mode
    
    input:
    path(fastqc_zips)  // all FastQC zip files collected
    
    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"),        emit: data
    
    script:
    """
    multiqc \\
        --force \\
        --filename multiqc_report.html \\
        --title "HiFi Raw Reads QC" \\
        --comment "Quality control metrics for HiFi reads" \\
        .
    """
    
    stub:
    """
    mkdir -p multiqc_report_data
    touch multiqc_report.html
    touch multiqc_report_data/multiqc_data.json
    """
}