/*
========================================================================================
    HI-C QC SUBWORKFLOW
========================================================================================
    Repo location: workflows/hic_qc.nf
    FastQC on Hi-C read pairs + MultiQC aggregation. Per-sample (meta).
========================================================================================
*/

nextflow.enable.dsl = 2

include { FASTQC_HIC } from '../modules/fastqc_hic'
include { MULTIQC_HIC } from '../modules/multiqc_hic'

workflow HIC_QC {
    take:
    hic_reads  // channel: tuple(meta, hic_r1, hic_r2)
    qc_label   // string: label for output directory (e.g. "raw" or "trimmed")

    main:
    // Run FastQC on each Hi-C read pair
    FASTQC_HIC(hic_reads, qc_label)

    // Aggregate all FastQC reports with MultiQC (zip files only)
    MULTIQC_HIC(
        FASTQC_HIC.out.fastqc_zip.map { meta, zips -> zips }.collect(),
        qc_label
    )

    emit:
    fastqc_html = FASTQC_HIC.out.fastqc_html
    fastqc_zip = FASTQC_HIC.out.fastqc_zip
    multiqc_report = MULTIQC_HIC.out.report
    multiqc_data = MULTIQC_HIC.out.data
}
