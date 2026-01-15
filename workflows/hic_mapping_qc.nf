/*
========================================================================================
    HI-C MAPPING QC SUBWORKFLOW
========================================================================================
    Comprehensive QC for Hi-C mapping to assemblies
    - Contact maps: visualization of Hi-C interactions
    - Pair statistics: valid pairs, trans/cis ratios, orientations
    - Coverage analysis: read distribution across assemblies
    - Combined report: aggregate all QC metrics
    
    Can process either raw or filtered BAMs with appropriate labeling
========================================================================================
*/

nextflow.enable.dsl = 2

include { HIC_CONTACT_MAP } from '../modules/hic_contact_map.nf'
include { HIC_PAIR_STATS } from '../modules/hic_pair_stats.nf'
include { HIC_COVERAGE } from '../modules/hic_coverage.nf'
include { COMBINE_HIC_MAPPING_QC } from '../modules/combine_hic_mapping_qc.nf'

workflow HIC_MAPPING_QC {
    take:
    bam_files      // channel: tuple(haplotype_id, bam, bai)
    assemblies     // channel: tuple(haplotype_id, assembly_fasta)
    qc_label       // string: "raw" or "filtered"
    
    main:
    
    /*
    ========================================================================================
        Add QC label to all channels
    ========================================================================================
    */
    bam_files
        .map { haplotype_id, bam, bai ->
            tuple(haplotype_id, bam, bai, qc_label)
        }
        .set { ch_bam_labeled }
    
    /*
    ========================================================================================
        Combine BAM files with assemblies for processes that need both
    ========================================================================================
    */
    ch_bam_labeled
        .map { haplotype_id, bam, bai, label ->
            tuple(haplotype_id, bam, bai, label)
        }
        .join(assemblies)
        .map { haplotype_id, bam, bai, label, assembly_fasta ->
            tuple(haplotype_id, bam, bai, assembly_fasta, label)
        }
        .set { ch_bam_with_assembly }
    
    /*
    ========================================================================================
        HI-C CONTACT MAP - Generate contact matrices and visualizations
    ========================================================================================
    */
    HIC_CONTACT_MAP(ch_bam_with_assembly)
    
    /*
    ========================================================================================
        HI-C PAIR STATISTICS - Analyze pair types, orientations, distances
    ========================================================================================
    */
    HIC_PAIR_STATS(ch_bam_labeled)
    
    /*
    ========================================================================================
        HI-C COVERAGE - Calculate and visualize coverage across assemblies
    ========================================================================================
    */
    HIC_COVERAGE(ch_bam_with_assembly)
    
    /*
    ========================================================================================
        COMBINE QC - Aggregate all Hi-C mapping QC results
    ========================================================================================
    */
    // Group all QC results by haplotype and qc_label
    
    HIC_CONTACT_MAP.out.stats
        .join(HIC_PAIR_STATS.out.summary, by: [0, 1])
        .join(HIC_COVERAGE.out.stats, by: [0, 1])
        .join(HIC_CONTACT_MAP.out.contact_maps.map { id, label, maps -> tuple(id, label, maps) }, by: [0, 1])
        .join(HIC_PAIR_STATS.out.plots.map { id, label, plots -> tuple(id, label, plots) }, by: [0, 1])
        .join(HIC_COVERAGE.out.plots.map { id, label, plots -> tuple(id, label, plots) }, by: [0, 1])
        .set { ch_all_hic_qc }
    
    COMBINE_HIC_MAPPING_QC(ch_all_hic_qc)
    
    emit:
    contact_maps = HIC_CONTACT_MAP.out.contact_maps
    cool_files = HIC_CONTACT_MAP.out.cool
    mcool_files = HIC_CONTACT_MAP.out.mcool
    pair_stats = HIC_PAIR_STATS.out.summary
    coverage_stats = HIC_COVERAGE.out.stats
    combined_report = COMBINE_HIC_MAPPING_QC.out.report
    combined_summary = COMBINE_HIC_MAPPING_QC.out.summary
}