/*
========================================================================================
    SUMMARY REPORT MODULE
========================================================================================
    Generates a comprehensive Markdown summary report for the genome assembly pipeline.
    
    Design:
    - ALL data passed via Nextflow channels (no publishDir scanning)
    - Markdown-only output (renders on GitHub/GitLab with proper image sizing)
    - Tables + collapsible tables for assembly QC metrics
    - Links to final genome assemblies, compiled QC CSV, and interactive HTML report
    - <img> tags with explicit width for snail plots, contact maps, dotplots
    - QC trend plots from COMPILE_FINAL_QC in collapsible section

    Publishes to ${params.outdir} directly so the report sits at the root of the
    results directory. All relative paths in the markdown are subdir/filename
    (e.g., snail_plots/sample_hap1_gap_filled_snail.svg).
========================================================================================
*/

process SUMMARY_REPORT {
    tag "summary_report"
    label 'summarize_assembly'

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path(report_manifest)       // TSV: type  id  id2  filename  subdir
    path(compiled_qc_csv)       // assembly_qc_metrics.csv from COMPILE_FINAL_QC
    path(telomere_summary)      // all_telomere_summaries.tsv  OR  NO_TELOMERES
    path(pairwise_summary)      // pairwise_alignment_summary.tsv  OR  NO_PAIRWISE

    output:
    path("assembly_report.md"), emit: report

    script:
    """
    Rscript ${projectDir}/r_scripts/generate_summary_report.R \
        --manifest ${report_manifest} \
        --compiled_qc ${compiled_qc_csv} \
        --telomere_summary ${telomere_summary} \
        --pairwise_summary ${pairwise_summary} \
        --output assembly_report.md
    """

    stub:
    """
    touch assembly_report.md
    """
}