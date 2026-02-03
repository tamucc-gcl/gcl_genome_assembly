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
    meryl_db     // channel: tuple(sample_id, meryl_db)
    qc_label     // value/channel: label for output subfolder
    
    main:
    
    /*
    ========================================================================================
        QUAST - Run on both haplotypes per sample
    ========================================================================================
    */
    QUAST(assemblies)
    
    /*
    ========================================================================================
        MERQURY - Run on both haplotypes per sample with pre-built meryl database
    ========================================================================================
    */
    assemblies
        .join(meryl_db)
        .map { sample_id, hap1_fasta, hap2_fasta, meryl_db ->
            tuple(sample_id, hap1_fasta, hap2_fasta, meryl_db)
        }
        .set { ch_merqury_input }
    
    MERQURY(ch_merqury_input)
    
    /*
    ========================================================================================
        Split haplotypes for per-haplotype QC
    ========================================================================================
    */
    assemblies
        .flatMap { sample_id, hap1_fasta, hap2_fasta ->
            [
                tuple("${sample_id}_hap1", sample_id, hap1_fasta),
                tuple("${sample_id}_hap2", sample_id, hap2_fasta)
            ]
        }
        .set { ch_individual_haplotypes }
    
    /*
    ========================================================================================
        BUSCO - Run on each haplotype independently
    ========================================================================================
    */
    BUSCO(
        ch_individual_haplotypes.map { haplotype_id, sample_id, fasta ->
            tuple(haplotype_id, fasta)
        }
    )
    
    /*
    ========================================================================================
        MAPPING QC - Map HiFi reads back to each haplotype
    ========================================================================================
    */
    ch_individual_haplotypes
        .map { haplotype_id, sample_id, fasta ->
            tuple(sample_id, haplotype_id, fasta)
        }
        .combine(hifi_reads, by: 0)
        .map { sample_id, haplotype_id, fasta, hifi_fastq ->
            tuple(haplotype_id, fasta, hifi_fastq)
        }
        .set { ch_mapping_input }
    
    MAPPING_QC(ch_mapping_input)
    
    /*
    ========================================================================================
        COMBINE QC - Aggregate all results and create plots
        
        CRITICAL FIX: Use groupTuple(by: 0) instead of groupTuple()
        - groupTuple() groups ALL samples together (changes when sample count changes)
        - groupTuple(by: 0) groups only by first element (sample_id) independently
    ========================================================================================
    */
    
    // Collect all QC outputs by sample
    QUAST.out.results
        .map { sample_id, results -> tuple(sample_id, results) }
        .set { ch_quast_by_sample }
    
    MERQURY.out.results
        .map { sample_id, results -> tuple(sample_id, results) }
        .set { ch_merqury_by_sample }
    
    // BUSCO results - FIX: groupTuple(by: 0) instead of groupTuple()
    BUSCO.out.results
        .map { haplotype_id, results -> 
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1] as Integer
            tuple(sample_id, hap_num, haplotype_id, results)
        }
        .groupTuple(by: 0)  // ← CRITICAL: Group only by sample_id (index 0)
        .map { sample_id, hap_nums, haplotype_ids, results_list ->
            def pairs = [hap_nums, haplotype_ids, results_list].transpose().sort { it[0] }
            tuple(sample_id, pairs.collect{it[1]}, pairs.collect{it[2]})
        }
        .set { ch_busco_by_sample }
    
    // MAPPING results - FIX: groupTuple(by: 0) instead of groupTuple()
    MAPPING_QC.out.results
        .map { haplotype_id, results ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1] as Integer
            tuple(sample_id, hap_num, haplotype_id, results)
        }
        .groupTuple(by: 0)  // ← CRITICAL: Group only by sample_id (index 0)
        .map { sample_id, hap_nums, haplotype_ids, results_list ->
            def pairs = [hap_nums, haplotype_ids, results_list].transpose().sort { it[0] }
            tuple(sample_id, pairs.collect{it[1]}, pairs.collect{it[2]})
        }
        .set { ch_mapping_by_sample }
    
    // Combine all QC results per sample
    ch_quast_by_sample
        .join(ch_merqury_by_sample)
        .join(ch_busco_by_sample)
        .join(ch_mapping_by_sample)
        .set { ch_all_qc }
    
    // Attach qc_label to each sample
    assemblies
        .map { sample_id, hap1_fasta, hap2_fasta -> tuple(sample_id, qc_label) }
        .set { ch_qc_label_by_sample }

    ch_all_qc
        .join(ch_qc_label_by_sample)
        .map { sample_id, quast_results, merqury_results, 
               haplotype_ids_busco, busco_results, 
               haplotype_ids_mapping, mapping_results, qc_label ->
            tuple(sample_id, qc_label,
                  quast_results, merqury_results,
                  haplotype_ids_busco, busco_results,
                  haplotype_ids_mapping, mapping_results)
        }
        .set { ch_all_qc_labeled }
        
    COMBINE_ASSEMBLY_QC(ch_all_qc_labeled)
    
    emit:
    assembly_summary = COMBINE_ASSEMBLY_QC.out.summary
}