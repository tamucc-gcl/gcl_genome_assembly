/*
========================================================================================
    HI-C SCAFFOLD QC WORKFLOW
========================================================================================
    Modular workflow for Hi-C QC on scaffolded assemblies via coordinate liftover
    
    Purpose:
    - Takes contig-coordinate pairs.gz + scaffold AGP + scaffold FASTA
    - Lifts pairs from contig to scaffold coordinates (NO REMAPPING)
    - Runs complete Hi-C QC in scaffold coordinate space
    
    This is a reusable "function" that can be called on:
    - Freshly scaffolded assemblies
    - Iterative scaffolding rounds
    - Manual scaffold edits with updated AGP
    
    Key Advantages:
    - No remapping needed - uses existing pairs.gz
    - Fast - AGP-based coordinate transformation
    - Accurate - preserves exact mapping positions from contigs
    - Can run multiple times with different scaffolds/AGPs
    
    Workflow Steps:
    1. Liftover: contig pairs.gz → scaffold pairs.gz (using AGP)
    2. QC: Run full Hi-C analysis in scaffold coordinates
    
    Outputs:
    - Lifted pairs.gz (scaffold coordinates)
    - Contact maps in scaffold coordinates
    - Pair statistics in scaffold coordinates
    - Coverage analysis in scaffold coordinates
    - Combined QC report
========================================================================================
*/

nextflow.enable.dsl = 2

include { HIC_LIFTOVER_PAIRS } from '../modules/hic_liftover_pairs.nf'
include { HIC_CONTACT_MAP_FROM_PAIRS } from '../modules/hic_contact_map.nf'
include { HIC_PAIR_STATS_FROM_PAIRS } from '../modules/hic_pair_stats.nf'
include { HIC_COVERAGE_FROM_PAIRS } from '../modules/hic_coverage.nf'
include { COMBINE_HIC_MAPPING_QC } from '../modules/combine_hic_mapping_qc.nf'

workflow HIC_SCAFFOLD_QC {
    take:
    contig_pairs       // channel: tuple(haplotype_id, pairs_gz) - in contig coordinates
    scaffold_agp       // channel: tuple(haplotype_id, agp) - contig→scaffold mapping
    scaffold_assemblies // channel: tuple(haplotype_id, scaffold_fasta)
    base_qc_label      // string: original label (e.g., "filtered")
    
    main:
    
    /*
    ========================================================================================
        Prepare inputs for liftover
        - Combine contig pairs with AGP and scaffold assembly
        - Prefix label with "scaffold_" to distinguish from contig QC
    ========================================================================================
    */
    contig_pairs
        .map { haplotype_id, label, pairs_gz ->
            // Remove any existing "scaffold_" prefix to avoid duplication
            def clean_label = label.replaceFirst(/^scaffold_/, '')
            tuple(haplotype_id, clean_label, pairs_gz)
        }
        .join(scaffold_agp, by: [0])
        .join(scaffold_assemblies, by: [0])
        .map { haplotype_id, label, pairs_gz, agp, scaffold_fasta ->
            // Add "scaffold_" prefix to distinguish scaffold QC outputs
            tuple(haplotype_id, "scaffold_${label}", pairs_gz, agp, scaffold_fasta)
        }
        .set { ch_liftover_input }
    
    /*
    // Debug: View liftover inputs
    ch_liftover_input.view { haplotype_id, label, pairs_gz, agp, scaffold_fasta ->
        """
        ========================================
        SCAFFOLD LIFTOVER INPUT
        Haplotype ID  : ${haplotype_id}
        Label         : ${label}
        Pairs (contig): ${pairs_gz}
        AGP           : ${agp}
        Scaffold FASTA: ${scaffold_fasta}
        ========================================
        """
    }
    */
    
    /*
    ========================================================================================
        STEP 1: Liftover pairs from contig to scaffold coordinates
    ========================================================================================
    */
    HIC_LIFTOVER_PAIRS(ch_liftover_input)
    
    /*
    ========================================================================================
        STEP 2: Prepare channels for QC in scaffold coordinates
    ========================================================================================
    */
    
    // For contact maps: need pairs + assembly + label
    HIC_LIFTOVER_PAIRS.out.pairs
        .join(scaffold_assemblies, by: [0])
        .map { haplotype_id, label, pairs_gz, scaffold_fasta ->
            tuple(haplotype_id, pairs_gz, scaffold_fasta, label)
        }
        .set { ch_scaffold_pairs_with_assembly }
    
    // For pair stats: need pairs + label
    HIC_LIFTOVER_PAIRS.out.pairs
        .map { haplotype_id, label, pairs_gz ->
            tuple(haplotype_id, pairs_gz, label)
        }
        .set { ch_scaffold_pairs_labeled }
    
    // For coverage: need pairs + assembly + label (same as contact maps)
    ch_scaffold_pairs_with_assembly
        .set { ch_scaffold_pairs_for_coverage }
    
    /*
    ========================================================================================
        STEP 3: Run complete Hi-C QC in scaffold coordinates
    ========================================================================================
    */
    
    // Contact maps
    HIC_CONTACT_MAP_FROM_PAIRS(ch_scaffold_pairs_with_assembly)
    
    // Pair statistics
    HIC_PAIR_STATS_FROM_PAIRS(ch_scaffold_pairs_labeled)
    
    // Coverage analysis
    HIC_COVERAGE_FROM_PAIRS(ch_scaffold_pairs_for_coverage)
    
    /*
    ========================================================================================
        STEP 4: Combine all QC results
    ========================================================================================
    */
    HIC_CONTACT_MAP_FROM_PAIRS.out.stats
        .join(HIC_PAIR_STATS_FROM_PAIRS.out.summary, by: [0, 1])
        .join(HIC_COVERAGE_FROM_PAIRS.out.stats, by: [0, 1])
        .join(HIC_CONTACT_MAP_FROM_PAIRS.out.contact_maps.map { id, label, maps -> tuple(id, label, maps) }, by: [0, 1])
        .join(HIC_PAIR_STATS_FROM_PAIRS.out.plots.map { id, label, plots -> tuple(id, label, plots) }, by: [0, 1])
        .join(HIC_COVERAGE_FROM_PAIRS.out.plots.map { id, label, plots -> tuple(id, label, plots) }, by: [0, 1])
        .set { ch_all_scaffold_qc }
    
    COMBINE_HIC_MAPPING_QC(ch_all_scaffold_qc)
    
    emit:
    // Lifted pairs for potential downstream use
    lifted_pairs = HIC_LIFTOVER_PAIRS.out.pairs
    
    // Core QC outputs (all in scaffold coordinates)
    contact_maps = HIC_CONTACT_MAP_FROM_PAIRS.out.contact_maps
    cool_files = HIC_CONTACT_MAP_FROM_PAIRS.out.cool
    mcool_files = HIC_CONTACT_MAP_FROM_PAIRS.out.mcool
    pair_stats = HIC_PAIR_STATS_FROM_PAIRS.out.summary
    coverage_stats = HIC_COVERAGE_FROM_PAIRS.out.stats
    combined_report = COMBINE_HIC_MAPPING_QC.out.report
    combined_summary = COMBINE_HIC_MAPPING_QC.out.summary
}