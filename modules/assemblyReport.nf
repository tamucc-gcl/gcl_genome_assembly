/*
========================================================================================
    ASSEMBLY REPORT MODULE
========================================================================================
    Generate an interactive HTML report comparing assembly QC metrics across
    all samples and pipeline stages. Automatically ranks assemblies by
    annotation suitability based on contiguity, base accuracy, gene space
    completeness, and scaffolding quality.

    The report is a fully self-contained HTML file (no external dependencies)
    with interactive filtering by stage, haplotype toggling, and visual
    metric comparisons.

    Inputs:
        - metrics_csv:  Collected QC metrics in long-format CSV
                        Columns: sample_id, stage, analysis, metric, hap1, hap2, both

    Outputs:
        - report_html:  Self-contained interactive HTML report

    Tools: Python 3 (stdlib only)
========================================================================================
*/

process ASSEMBLY_REPORT {
    tag "report"
    label 'assembly_report'
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

    input:
    path(metrics_csv)

    output:
    path("assembly_qc_report.html"), emit: report_html

    script:
    """
    python3 ${projectDir}/py_scripts/generate_assembly_report.py \\
        --input ${metrics_csv} \\
        --output assembly_qc_report.html \\
        --report-stage "${params.report_stage ?: 'final'}"
    """

    stub:
    """
    touch assembly_qc_report.html
    """
}