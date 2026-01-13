/*
========================================================================================
    MULTIQC HI-C MODULE
========================================================================================
    Aggregates FastQC reports from Hi-C reads into a single MultiQC report
========================================================================================
*/

process MULTIQC_HIC {
    label 'multiqc'
    
    publishDir "${params.outdir}/qc/hic", mode: params.publish_dir_mode
    
    input:
    path(fastqc_zips)  // all FastQC zip files collected
    
    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"),        emit: data
    
    script:
    """
    multiqc \\
        --force \\
        --filename multiqc_report.html \\
        --title "Hi-C Raw Reads QC" \\
        --comment "Quality control metrics for Hi-C paired-end reads" \\
        .
    """
    
    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    touch multiqc_data/multiqc_data.json
    """
}