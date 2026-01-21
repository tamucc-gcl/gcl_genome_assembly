/*
========================================================================================
    HI-C QC FROM PAIRS WORKFLOW
========================================================================================
    Modular workflow for Hi-C QC starting from pairs.gz files
    
    Purpose:
    - Takes pre-computed pairs.gz files and assemblies as input
    - Runs complete Hi-C QC analysis WITHOUT remapping
    - Useful for QC on lifted coordinates, gap-filled assemblies, etc.
    
    This is a reusable "function" that can be called on:
    - Lifted pairs (scaffold coordinates)
    - Pairs from external sources
    - Reprocessed pairs with different filtering
    - Gap-filled assembly coordinates
    
    Key Advantage:
    - No remapping needed - much faster for iterative assembly improvement
    - Can run QC on multiple assembly versions using the same pairs
    
    Outputs:
    - Contact maps and visualizations
    - Pair statistics (cis/trans ratios, orientations, etc.)
    - Coverage analysis
    - Combined QC report
========================================================================================
*/

nextflow.enable.dsl = 2

include { HIC_CONTACT_MAP_FROM_PAIRS } from '../modules/hic_contact_map.nf'
include { HIC_PAIR_STATS_FROM_PAIRS } from '../modules/hic_pair_stats.nf'
include { HIC_COVERAGE_FROM_PAIRS } from '../modules/hic_coverage.nf'
include { COMBINE_HIC_MAPPING_QC } from '../modules/combine_hic_mapping_qc.nf'

workflow HIC_QC_FROM_PAIRS {
    take:
    pairs_files    // channel: tuple(haplotype_id, pairs_gz)
    assemblies     // channel: tuple(haplotype_id, assembly_fasta)
    qc_label       // string: "scaffold_filtered", "gapfilled", or custom label
    
    main:
    
    /*
    ========================================================================================
        Add QC label to pairs files
    ========================================================================================
    */
    pairs_files
        .map { haplotype_id, pairs_gz ->
            tuple(haplotype_id, pairs_gz, qc_label)
        }
        .set { ch_pairs_labeled }
    
    /*
    ========================================================================================
        Combine pairs files with assemblies
    ========================================================================================
    */
    pairs_files
        .join(assemblies)
        .map { haplotype_id, pairs_gz, assembly_fasta ->
            tuple(haplotype_id, pairs_gz, assembly_fasta, qc_label)
        }
        .set { ch_pairs_with_assembly }
    
    /*
    ========================================================================================
        HI-C CONTACT MAP - Generate contact matrices from pairs.gz
    ========================================================================================
    */
    HIC_CONTACT_MAP_FROM_PAIRS(ch_pairs_with_assembly)
    
    /*
    ========================================================================================
        HI-C PAIR STATISTICS - Analyze pair types from pairs.gz
    ========================================================================================
    */
    HIC_PAIR_STATS_FROM_PAIRS(ch_pairs_labeled)
    
    /*
    ========================================================================================
        HI-C COVERAGE - Calculate coverage from pairs.gz
    ========================================================================================
    */
    HIC_COVERAGE_FROM_PAIRS(ch_pairs_with_assembly)
    
    /*
    ========================================================================================
        COMBINE QC - Aggregate all Hi-C QC results
    ========================================================================================
    */
    HIC_CONTACT_MAP_FROM_PAIRS.out.stats
        .join(HIC_PAIR_STATS_FROM_PAIRS.out.summary, by: [0, 1])
        .join(HIC_COVERAGE_FROM_PAIRS.out.stats, by: [0, 1])
        .join(HIC_CONTACT_MAP_FROM_PAIRS.out.contact_maps.map { id, label, maps -> tuple(id, label, maps) }, by: [0, 1])
        .join(HIC_PAIR_STATS_FROM_PAIRS.out.plots.map { id, label, plots -> tuple(id, label, plots) }, by: [0, 1])
        .join(HIC_COVERAGE_FROM_PAIRS.out.plots.map { id, label, plots -> tuple(id, label, plots) }, by: [0, 1])
        .set { ch_all_qc }
    
    COMBINE_HIC_MAPPING_QC(ch_all_qc)
    
    emit:
    // Core QC outputs
    contact_maps = HIC_CONTACT_MAP_FROM_PAIRS.out.contact_maps
    cool_files = HIC_CONTACT_MAP_FROM_PAIRS.out.cool
    mcool_files = HIC_CONTACT_MAP_FROM_PAIRS.out.mcool
    pair_stats = HIC_PAIR_STATS_FROM_PAIRS.out.summary
    coverage_stats = HIC_COVERAGE_FROM_PAIRS.out.stats
    combined_report = COMBINE_HIC_MAPPING_QC.out.report
    combined_summary = COMBINE_HIC_MAPPING_QC.out.summary
}