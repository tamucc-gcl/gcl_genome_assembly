/*
========================================================================================
    HiFi QC SUBWORKFLOW
========================================================================================
    Runs FastQC on HiFi reads and aggregates results with MultiQC
========================================================================================
*/

nextflow.enable.dsl = 2

include { FASTQC_HIFI } from '../modules/fastqc_hifi.nf'
include { MULTIQC_HIFI } from '../modules/multiqc_hifi.nf'

workflow HIFI_QC {
    take:
    hifi_reads  // channel: tuple(sample_id, hic_r1, hic_r2)
    
    main:
    // Run FastQC on each Hi-C read pair
    FASTQC_HIFI(hifi_reads)
    
    // Aggregate all FastQC reports with MultiQC
    // Extract just the zip files from the tuples
    MULTIQC_HIFI(
        FASTQC_HIFI.out.fastqc_zip.map { sample_id, zips -> zips }.collect()
    )
    
    emit:
    fastqc_html = FASTQC_HIFI.out.fastqc_html
    fastqc_zip = FASTQC_HIFI.out.fastqc_zip
    multiqc_report = MULTIQC_HIFI.out.report
    multiqc_data = MULTIQC_HIFI.out.data
}