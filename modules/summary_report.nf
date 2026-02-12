/*
========================================================================================
    SUMMARY REPORT MODULE
========================================================================================
    Generates a markdown summary report with:
    1. Visual table with links to snail plots, contact maps, and dotplots
    2. QC metrics table from compiled TSVs
    3. Links to all output files
    
    Relies on files already published to params.outdir by other processes
========================================================================================
*/

process SUMMARY_REPORT {
    tag "summary_report"
    label 'summarize_assembly'
    
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode
    
    input:
    path(qc_inputs_dir)  // From COMPILE_FINAL_QC.out.inputs_dir
    
    output:
    path("pipeline_summary_report.md"), emit: md_report
    path("pipeline_summary_report.html"), emit: html_report
    
    script:
    """
    Rscript ${projectDir}/r_scripts/generate_summary_report.R \\
        --qc_dir ${qc_inputs_dir} \\
        --outdir_base "${params.outdir}" \\
        --output_dir .
    """
    
    stub:
    """
    touch pipeline_summary_report.md
    touch pipeline_summary_report.html
    """
}