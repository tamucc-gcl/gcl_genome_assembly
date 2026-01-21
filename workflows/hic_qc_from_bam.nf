/*
========================================================================================
    HI-C QC FROM BAM WORKFLOW
========================================================================================
    Modular workflow for Hi-C QC starting from BAM files
    
    Purpose:
    - Takes BAM files and assemblies as input
    - Runs complete Hi-C QC analysis
    - Outputs pairs.gz for potential downstream use (e.g., scaffold liftover)
    
    This is a reusable "function" that can be called on:
    - Raw mapped BAMs (qc_label = "raw")
    - Filtered BAMs (qc_label = "filtered")
    - Any other BAM with custom label
    
    Outputs:
    - Contact maps and visualizations
    - Pair statistics (cis/trans ratios, orientations, etc.)
    - Coverage analysis
    - Combined QC report
    - pairs.gz file (for downstream liftover or other analyses)
========================================================================================
*/

nextflow.enable.dsl = 2

include { HIC_CONTACT_MAP } from '../modules/hic_contact_map.nf'
include { HIC_PAIR_STATS } from '../modules/hic_pair_stats.nf'
include { HIC_COVERAGE } from '../modules/hic_coverage.nf'
include { COMBINE_HIC_MAPPING_QC } from '../modules/combine_hic_mapping_qc.nf'

workflow HIC_QC_FROM_BAM {
    take:
    bam_files      // channel: tuple(haplotype_id, bam, bai)
    assemblies     // channel: tuple(haplotype_id, assembly_fasta)
    qc_label       // string: "raw", "filtered", or custom label
    
    main:
    
    /*
    ========================================================================================
        Add QC label to BAM files
    ========================================================================================
    */
    bam_files
        .map { haplotype_id, bam, bai ->
            tuple(haplotype_id, bam, bai, qc_label)
        }
        .set { ch_bam_labeled }
    
    /*
    ========================================================================================
        Combine BAM files with assemblies
    ========================================================================================
    */
    ch_bam_labeled
        .join(assemblies)
        .map { haplotype_id, bam, bai, label, assembly_fasta ->
            tuple(haplotype_id, bam, bai, assembly_fasta, label)
        }
        .set { ch_bam_with_assembly }
    
    /*
    ========================================================================================
        HI-C CONTACT MAP - Generate contact matrices and visualizations
        Also produces pairs.gz for downstream use
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
        COMBINE QC - Aggregate all Hi-C QC results
    ========================================================================================
    */
    HIC_CONTACT_MAP.out.stats
        .join(HIC_PAIR_STATS.out.summary, by: [0, 1])
        .join(HIC_COVERAGE.out.stats, by: [0, 1])
        .join(HIC_CONTACT_MAP.out.contact_maps.map { id, label, maps -> tuple(id, label, maps) }, by: [0, 1])
        .join(HIC_PAIR_STATS.out.plots.map { id, label, plots -> tuple(id, label, plots) }, by: [0, 1])
        .join(HIC_COVERAGE.out.plots.map { id, label, plots -> tuple(id, label, plots) }, by: [0, 1])
        .set { ch_all_qc }
    
    COMBINE_HIC_MAPPING_QC(ch_all_qc)
    
    emit:
    // Core QC outputs
    contact_maps = HIC_CONTACT_MAP.out.contact_maps
    cool_files = HIC_CONTACT_MAP.out.cool
    mcool_files = HIC_CONTACT_MAP.out.mcool
    pair_stats = HIC_PAIR_STATS.out.summary
    coverage_stats = HIC_COVERAGE.out.stats
    combined_report = COMBINE_HIC_MAPPING_QC.out.report
    combined_summary = COMBINE_HIC_MAPPING_QC.out.summary
    
    // pairs.gz for downstream use (e.g., liftover to scaffold coordinates)
    pairs = HIC_CONTACT_MAP.out.pairs
}