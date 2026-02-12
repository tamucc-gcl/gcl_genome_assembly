/*
========================================================================================
    SUMMARY REPORT MODULE
========================================================================================
    Generates a markdown summary report by scanning published output directories.
    Takes no inputs from other processes - just needs params.outdir.
========================================================================================
*/

process SUMMARY_REPORT {
    tag "summary_report"
    label 'summarize_assembly'
    
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode
    
    input:
    val(trigger)  // Just a signal that upstream processes are done
    
    output:
    path("pipeline_summary_report.md"), emit: md_report
    path("pipeline_summary_report.html"), emit: html_report
    
    script:
    def abs_outdir = file(params.outdir).toAbsolutePath()
    """
    Rscript ${projectDir}/r_scripts/generate_summary_report.R \\
        --outdir "${abs_outdir}" \\
        --output_dir .
    """
    
    stub:
    """
    touch pipeline_summary_report.md
    touch pipeline_summary_report.html
    """
}