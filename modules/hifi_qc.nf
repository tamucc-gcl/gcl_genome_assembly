/*
========================================================================================
    HiFi QC SUBWORKFLOW
========================================================================================
    Repo location: workflows/hifi_qc.nf
    FastQC on HiFi reads + MultiQC aggregation. Per-sample (meta).
========================================================================================
*/

nextflow.enable.dsl = 2

include { FASTQC_HIFI } from '../modules/fastqc_hifi.nf'
include { MULTIQC_HIFI } from '../modules/multiqc_hifi.nf'

workflow HIFI_QC {
    take:
    hifi_reads  // channel: tuple(meta, hifi_fastq)

    main:
    // Run FastQC on each HiFi read set
    FASTQC_HIFI(hifi_reads)

    // Aggregate all FastQC reports with MultiQC (zip files only)
    MULTIQC_HIFI(
        FASTQC_HIFI.out.fastqc_zip.map { meta, zips -> zips }.collect()
    )

    emit:
    fastqc_html = FASTQC_HIFI.out.fastqc_html
    fastqc_zip = FASTQC_HIFI.out.fastqc_zip
    multiqc_report = MULTIQC_HIFI.out.report
    multiqc_data = MULTIQC_HIFI.out.data
}
