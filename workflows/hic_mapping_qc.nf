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

include { HIC_CONTACT_MAP; HIC_CONTACT_MAP_FROM_PAIRS } from '../modules/hic_contact_map.nf'
include { HIC_PAIR_STATS; HIC_PAIR_STATS_FROM_PAIRS } from '../modules/hic_pair_stats.nf'
include { HIC_COVERAGE; HIC_COVERAGE_FROM_PAIRS } from '../modules/hic_coverage.nf'
include { COMBINE_HIC_MAPPING_QC as COMBINE_HIC_MAPPING_QC_CONTIG } from '../modules/combine_hic_mapping_qc.nf'
include { COMBINE_HIC_MAPPING_QC as COMBINE_HIC_MAPPING_QC_SCAFFOLD } from '../modules/combine_hic_mapping_qc.nf'
include { HIC_LIFTOVER_PAIRS } from '../modules/hic_liftover_pairs.nf'

workflow HIC_MAPPING_QC {
    take:
    bam_files      // channel: tuple(haplotype_id, bam, bai)
    assemblies     // channel: tuple(haplotype_id, assembly_fasta)
    qc_label       // string: "raw" or "filtered"
    scaffold_assemblies   // channel: tuple(haplotype_id, scaffold_fasta)
    scaffold_agp          // channel: tuple(haplotype_id, scaffold_agp)
    
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
    
    COMBINE_HIC_MAPPING_QC_CONTIG(ch_all_hic_qc)


/*
========================================================================================
    OPTIONAL: SCAFFOLD-COORDINATE QC (NO REMAPPING)
========================================================================================
    Uses pairs.gz from contig-mapped BAM (HIC_CONTACT_MAP output), lifts those pairs to
    scaffold coordinates using the scaffold AGP, then re-runs pairtools/cooler QC on the
    scaffold assembly without remapping.

    Requirements:
    - HIC_CONTACT_MAP must emit pairs.gz
    - scaffold_assemblies: FASTA of scaffolded assembly for each haplotype
    - scaffold_agp: AGP describing contig->scaffold placement for each haplotype
========================================================================================
*/

// 1) Lift contig pairs.gz to scaffold coordinates
HIC_CONTACT_MAP.out.pairs
    .join(scaffold_agp, by: [0])          // join on haplotype_id
    .join(scaffold_assemblies, by: [0])   // join on haplotype_id
    .map { haplotype_id, label, pairs_gz, agp, scaffold_fasta ->
        tuple(haplotype_id, "scaffold_${label}", pairs_gz, agp, scaffold_fasta)
    }
    .set { ch_pairs_for_liftover }

HIC_LIFTOVER_PAIRS(ch_pairs_for_liftover)

// 2) Run scaffold-coordinate QC starting from lifted pairs.gz
HIC_LIFTOVER_PAIRS.out.pairs
    .join(scaffold_assemblies, by: [0])
    .map { haplotype_id, label, pairs_gz, scaffold_fasta ->
        tuple(haplotype_id, pairs_gz, scaffold_fasta, label)
    }
    .set { ch_scaffold_pairs_with_assembly }

HIC_CONTACT_MAP_FROM_PAIRS(ch_scaffold_pairs_with_assembly)

HIC_LIFTOVER_PAIRS.out.pairs
    .map { haplotype_id, label, pairs_gz ->
        tuple(haplotype_id, pairs_gz, label)
    }
    .set { ch_scaffold_pairs_labeled }

HIC_PAIR_STATS_FROM_PAIRS(ch_scaffold_pairs_labeled)

HIC_LIFTOVER_PAIRS.out.pairs
    .join(scaffold_assemblies, by: [0])
    .map { haplotype_id, label, pairs_gz, scaffold_fasta ->
        tuple(haplotype_id, pairs_gz, scaffold_fasta, label)
    }
    .set { ch_scaffold_pairs_for_coverage }

HIC_COVERAGE_FROM_PAIRS(ch_scaffold_pairs_for_coverage)

/*
========================================================================================
    COMBINE QC (SCAFFOLD COORDS)
========================================================================================
*/
HIC_CONTACT_MAP_FROM_PAIRS.out.stats
    .join(HIC_PAIR_STATS_FROM_PAIRS.out.summary, by: [0, 1])
    .join(HIC_COVERAGE_FROM_PAIRS.out.stats, by: [0, 1])
    .join(HIC_CONTACT_MAP_FROM_PAIRS.out.contact_maps.map { id, label, maps -> tuple(id, label, maps) }, by: [0, 1])
    .join(HIC_PAIR_STATS_FROM_PAIRS.out.plots.map { id, label, plots -> tuple(id, label, plots) }, by: [0, 1])
    .join(HIC_COVERAGE_FROM_PAIRS.out.plots.map { id, label, plots -> tuple(id, label, plots) }, by: [0, 1])
    .set { ch_all_hic_qc_scaffold }

COMBINE_HIC_MAPPING_QC_SCAFFOLD(ch_all_hic_qc_scaffold)
    
    emit:
    contact_maps = HIC_CONTACT_MAP.out.contact_maps
    cool_files = HIC_CONTACT_MAP.out.cool
    mcool_files = HIC_CONTACT_MAP.out.mcool
    pairs = HIC_CONTACT_MAP.out.pairs
    pair_stats = HIC_PAIR_STATS.out.summary
    coverage_stats = HIC_COVERAGE.out.stats
    combined_report = COMBINE_HIC_MAPPING_QC_CONTIG.out.report
    combined_summary = COMBINE_HIC_MAPPING_QC_CONTIG.out.summary
    // Scaffold-coordinate QC outputs (labels prefixed with "scaffold_")
    scaffold_contact_maps = HIC_CONTACT_MAP_FROM_PAIRS.out.contact_maps
    scaffold_cool_files = HIC_CONTACT_MAP_FROM_PAIRS.out.cool
    scaffold_mcool_files = HIC_CONTACT_MAP_FROM_PAIRS.out.mcool
    scaffold_pair_stats = HIC_PAIR_STATS_FROM_PAIRS.out.summary
    scaffold_coverage_stats = HIC_COVERAGE_FROM_PAIRS.out.stats
    scaffold_combined_report = COMBINE_HIC_MAPPING_QC_SCAFFOLD.out.report
    scaffold_combined_summary = COMBINE_HIC_MAPPING_QC_SCAFFOLD.out.summary
}