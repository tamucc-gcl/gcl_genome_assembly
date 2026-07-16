/*
========================================================================================
    SHORT-READ QC SUBWORKFLOW
========================================================================================
    Repo location: workflows/shortread_qc.nf

    FastQC on short-read (Illumina) paired-end reads + MultiQC aggregation. Per-sample
    (meta). Mirrors HIC_QC — the same input-read QC tools applied to the short-read inputs.
    FASTQC_HIC is read-type-agnostic (runs FastQC on two FASTQs), so it's reused via alias;
    MULTIQC needs short-read titles/path, hence a dedicated MULTIQC_SHORTREAD.
========================================================================================
*/

nextflow.enable.dsl = 2

include { FASTQC_HIC as FASTQC_SHORTREAD } from '../modules/fastqc_hic'
include { MULTIQC_SHORTREAD }              from '../modules/multiqc_shortread'

workflow SHORTREAD_QC {
    take:
    sr_reads   // channel: tuple(meta, sr_r1, sr_r2)
    qc_label   // string: label for output directory (e.g. "raw")

    main:
    // FastQC on each short-read pair
    FASTQC_SHORTREAD(sr_reads, qc_label)

    // Aggregate all FastQC reports with MultiQC (zip files only)
    MULTIQC_SHORTREAD(
        FASTQC_SHORTREAD.out.fastqc_zip.map { meta, zips -> zips }.collect(),
        qc_label
    )

    emit:
    fastqc_html    = FASTQC_SHORTREAD.out.fastqc_html
    fastqc_zip     = FASTQC_SHORTREAD.out.fastqc_zip
    multiqc_report = MULTIQC_SHORTREAD.out.report
    multiqc_data   = MULTIQC_SHORTREAD.out.data
}
