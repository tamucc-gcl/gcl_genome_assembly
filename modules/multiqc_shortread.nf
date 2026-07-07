/*
========================================================================================
    MULTIQC SHORT-READ MODULE
========================================================================================
    Repo location: modules/multiqc_shortread.nf

    Aggregates FastQC reports from short-read (Illumina) paired-end reads into a single
    MultiQC report. Mirrors MULTIQC_HIC but with short-read titles and its own output path.
========================================================================================
*/

process MULTIQC_SHORTREAD {
    label 'multiqc'

    publishDir "${params.outdir}/qc/shortread/${qc_label}", mode: params.publish_dir_mode

    input:
    path(fastqc_zips)  // all FastQC zip files collected
    val qc_label       // e.g. "raw"

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data

    script:
    def title = qc_label == "raw" ? "Short-read Raw Reads QC" : "Short-read Trimmed Reads QC"
    """
    multiqc \\
        --force \\
        --filename multiqc_report.html \\
        --title "${title}" \\
        --comment "Quality control metrics for short-read (Illumina) paired-end reads (${qc_label})" \\
        .
    """

    stub:
    """
    mkdir -p multiqc_report_data
    touch multiqc_report.html
    touch multiqc_report_data/multiqc_data.json
    """
}
