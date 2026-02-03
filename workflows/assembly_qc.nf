/*
========================================================================================
    ASSEMBLY QC SUBWORKFLOW (CACHE-OPTIMIZED)
========================================================================================
    Comprehensive QC for phased genome assemblies
    - QUAST: per-sample (both haplotypes)
    - MERQURY: per-sample (both haplotypes) - uses pre-built meryl database
    - BUSCO: per-haplotype
    - MAPPING: per-haplotype (map HiFi reads back to assembly)
    - COMBINE_QC: aggregate and visualize all QC metrics
    
    CACHING FIX:
    - Each sample processes completely independently
    - No cross-sample channel operations that break cache
    - Adding new samples won't invalidate existing sample caches
========================================================================================
*/

nextflow.enable.dsl = 2

include { QUAST } from '../modules/quast.nf'
include { MERQURY } from '../modules/merqury.nf'
include { BUSCO } from '../modules/busco.nf'
include { MAPPING_QC } from '../modules/mapping_qc.nf'
include { COMBINE_ASSEMBLY_QC } from '../modules/combine_assembly_qc.nf'

workflow ASSEMBLY_QC {
    take:
    assemblies   // channel: tuple(sample_id, hap1_fasta, hap2_fasta)
    hifi_reads   // channel: tuple(sample_id, hifi_fastq)
    meryl_db     // channel: tuple(sample_id, meryl_db) - PRE-BUILT k-mer database
    qc_label     // value: label for output subfolder (e.g. 'contig' or 'scaffold')
    
    main:
    
    /*
    ========================================================================================
        CRITICAL: Add qc_label to assemblies channel FIRST
        This ensures each (sample_id, qc_label) combination is independently cacheable
    ========================================================================================
    */
    assemblies
        .map { sample_id, hap1_fasta, hap2_fasta ->
            tuple(sample_id, qc_label, hap1_fasta, hap2_fasta)
        }
        .set { ch_assemblies_labeled }
    
    /*
    ========================================================================================
        QUAST - Run on both haplotypes per sample
        Input: (sample_id, qc_label, hap1_fasta, hap2_fasta) - stable per sample
    ========================================================================================
    */
    ch_assemblies_labeled
        .map { sample_id, label, hap1_fasta, hap2_fasta ->
            tuple(sample_id, hap1_fasta, hap2_fasta)
        }
        .set { ch_quast_input }
    
    QUAST(ch_quast_input)
    
    /*
    ========================================================================================
        MERQURY - Run on both haplotypes per sample with pre-built meryl database
        Input: (sample_id, hap1_fasta, hap2_fasta, meryl_db) - stable per sample
    ========================================================================================
    */
    ch_assemblies_labeled
        .map { sample_id, label, hap1_fasta, hap2_fasta ->
            tuple(sample_id, hap1_fasta, hap2_fasta)
        }
        .join(meryl_db, by: 0)
        .map { sample_id, hap1_fasta, hap2_fasta, meryl_db ->
            tuple(sample_id, hap1_fasta, hap2_fasta, meryl_db)
        }
        .set { ch_merqury_input }
    
    MERQURY(ch_merqury_input)
    
    /*
    ========================================================================================
        Split haplotypes for per-haplotype QC
        CRITICAL: Use flatMap on the labeled channel to maintain sample independence
    ========================================================================================
    */
    ch_assemblies_labeled
        .flatMap { sample_id, label, hap1_fasta, hap2_fasta ->
            [
                tuple(sample_id, "${sample_id}_hap1", hap1_fasta),
                tuple(sample_id, "${sample_id}_hap2", hap2_fasta)
            ]
        }
        .set { ch_individual_haplotypes }
    
    /*
    ========================================================================================
        BUSCO - Run on each haplotype independently
        Input: (haplotype_id, fasta) - stable per haplotype
    ========================================================================================
    */
    ch_individual_haplotypes
        .map { sample_id, haplotype_id, fasta ->
            tuple(haplotype_id, fasta)
        }
        .set { ch_busco_input }
    
    BUSCO(ch_busco_input)
    
    /*
    ========================================================================================
        MAPPING QC - Map HiFi reads back to each haplotype
        Input: (haplotype_id, fasta, hifi_fastq) - stable per haplotype
    ========================================================================================
    */
    ch_individual_haplotypes
        .combine(hifi_reads, by: 0)
        .map { sample_id, haplotype_id, fasta, hifi_fastq ->
            tuple(haplotype_id, fasta, hifi_fastq)
        }
        .set { ch_mapping_input }
    
    MAPPING_QC(ch_mapping_input)
    
    /*
    ========================================================================================
        COMBINE QC - Aggregate all results per sample
        CRITICAL FIX: Join operations based on sample_id, not groupTuple
        This ensures each sample's combination is stable and independently cacheable
    ========================================================================================
    */
    
    // Start with labeled assemblies to ensure we have (sample_id, qc_label) for each
    ch_assemblies_labeled
        .map { sample_id, label, hap1, hap2 -> tuple(sample_id, label) }
        .set { ch_sample_labels }
    
    // Join QUAST results
    ch_sample_labels
        .join(QUAST.out.results, by: 0)
        .set { ch_with_quast }
    
    // Join MERQURY results
    ch_with_quast
        .join(MERQURY.out.results, by: 0)
        .set { ch_with_merqury }
    
    // Collect BUSCO results per sample (must maintain stable order: hap1, hap2)
    BUSCO.out.results
        .map { haplotype_id, results ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1] as Integer
            tuple(sample_id, hap_num, haplotype_id, results)
        }
        .groupTuple(by: 0)  // Group by sample_id only
        .map { sample_id, hap_nums, haplotype_ids, results_list ->
            // Sort by haplotype number to ensure stable hap1, hap2 order
            def sorted = [hap_nums, haplotype_ids, results_list].transpose().sort { it[0] }
            tuple(sample_id, sorted.collect{it[1]}, sorted.collect{it[2]})
        }
        .set { ch_busco_by_sample }
    
    // Join BUSCO results
    ch_with_merqury
        .join(ch_busco_by_sample, by: 0)
        .set { ch_with_busco }
    
    // Collect MAPPING results per sample (must maintain stable order: hap1, hap2)
    MAPPING_QC.out.results
        .map { haplotype_id, results ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1] as Integer
            tuple(sample_id, hap_num, haplotype_id, results)
        }
        .groupTuple(by: 0)  // Group by sample_id only
        .map { sample_id, hap_nums, haplotype_ids, results_list ->
            // Sort by haplotype number to ensure stable hap1, hap2 order
            def sorted = [hap_nums, haplotype_ids, results_list].transpose().sort { it[0] }
            tuple(sample_id, sorted.collect{it[1]}, sorted.collect{it[2]})
        }
        .set { ch_mapping_by_sample }
    
    // Join MAPPING results
    ch_with_busco
        .join(ch_mapping_by_sample, by: 0)
        .map { sample_id, label, quast_results, merqury_results, 
               haplotype_ids_busco, busco_results,
               haplotype_ids_mapping, mapping_results ->
            tuple(sample_id, label, quast_results, merqury_results,
                  haplotype_ids_busco, busco_results,
                  haplotype_ids_mapping, mapping_results)
        }
        .set { ch_all_qc }
    
    COMBINE_ASSEMBLY_QC(ch_all_qc)
    
    emit:
    assembly_summary = COMBINE_ASSEMBLY_QC.out.summary
}