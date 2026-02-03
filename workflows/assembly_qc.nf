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
    qc_label     // value: label for output subfolder (e.g. 'contig' or 'scaffold')
    
    main:
    
    /*
    ========================================================================================
        CRITICAL: Tag each sample with qc_label IMMEDIATELY
        This makes (sample_id, qc_label) the stable cache key throughout
    ========================================================================================
    */
    assemblies
        .map { sample_id, hap1_fasta, hap2_fasta ->
            tuple(sample_id, qc_label, hap1_fasta, hap2_fasta)
        }
        .set { ch_assemblies_tagged }
    
    /*
    ========================================================================================
        QUAST - Run on both haplotypes per sample
    ========================================================================================
    */
    QUAST(
        ch_assemblies_tagged.map { sample_id, label, hap1, hap2 ->
            tuple(sample_id, hap1, hap2)
        }
    )
    
    /*
    ========================================================================================
        MERQURY - Run on both haplotypes per sample with pre-built meryl database
    ========================================================================================
    */
    ch_assemblies_tagged
        .map { sample_id, label, hap1_fasta, hap2_fasta ->
            tuple(sample_id, hap1_fasta, hap2_fasta)
        }
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
    ch_assemblies_tagged
        .flatMap { sample_id, label, hap1_fasta, hap2_fasta ->
            [
                tuple("${sample_id}_hap1", sample_id, label, hap1_fasta),
                tuple("${sample_id}_hap2", sample_id, label, hap2_fasta)
            ]
        }
        .set { ch_individual_haplotypes }
    
    /*
    ========================================================================================
        BUSCO - Run on each haplotype independently
    ========================================================================================
    */
    BUSCO(
        ch_individual_haplotypes.map { haplotype_id, sample_id, label, fasta ->
            tuple(haplotype_id, fasta)
        }
    )
    
    /*
    ========================================================================================
        MAPPING QC - Map HiFi reads back to each haplotype
    ========================================================================================
    */
    ch_individual_haplotypes
        .map { haplotype_id, sample_id, label, fasta ->
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
        COMBINE QC - Aggregate all results per (sample_id, qc_label)
        
        KEY INSIGHT: We need to preserve qc_label through all the grouping operations
        so that COMBINE_ASSEMBLY_QC gets (sample_id, qc_label, ...) as input
    ========================================================================================
    */
    
    // Tag QUAST results with qc_label
    ch_assemblies_tagged
        .map { sample_id, label, hap1, hap2 -> tuple(sample_id, label) }
        .join(QUAST.out.results)
        .map { sample_id, label, results -> tuple(sample_id, label, results) }
        .set { ch_quast_tagged }
    
    // Tag MERQURY results with qc_label  
    ch_assemblies_tagged
        .map { sample_id, label, hap1, hap2 -> tuple(sample_id, label) }
        .join(MERQURY.out.results)
        .map { sample_id, label, results -> tuple(sample_id, label, results) }
        .set { ch_merqury_tagged }
    
    // Collect and tag BUSCO results - groupTuple(by: 0) groups only by sample_id
    BUSCO.out.results
        .map { haplotype_id, results ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1] as Integer
            tuple(sample_id, hap_num, haplotype_id, results)
        }
        .groupTuple(by: 0, sort: true)  // Group only by sample_id
        .map { sample_id, hap_nums, haplotype_ids, results_list ->
            def pairs = [hap_nums, haplotype_ids, results_list].transpose().sort { it[0] }
            tuple(sample_id, pairs.collect{it[1]}, pairs.collect{it[2]})
        }
        .set { ch_busco_grouped }
    
    // Tag with qc_label
    ch_assemblies_tagged
        .map { sample_id, label, hap1, hap2 -> tuple(sample_id, label) }
        .join(ch_busco_grouped)
        .map { sample_id, label, haplotype_ids, results ->
            tuple(sample_id, label, haplotype_ids, results)
        }
        .set { ch_busco_tagged }
    
    // Collect and tag MAPPING results - groupTuple(by: 0) groups only by sample_id
    MAPPING_QC.out.results
        .map { haplotype_id, results ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1] as Integer
            tuple(sample_id, hap_num, haplotype_id, results)
        }
        .groupTuple(by: 0, sort: true)  // Group only by sample_id
        .map { sample_id, hap_nums, haplotype_ids, results_list ->
            def pairs = [hap_nums, haplotype_ids, results_list].transpose().sort { it[0] }
            tuple(sample_id, pairs.collect{it[1]}, pairs.collect{it[2]})
        }
        .set { ch_mapping_grouped }
    
    // Tag with qc_label
    ch_assemblies_tagged
        .map { sample_id, label, hap1, hap2 -> tuple(sample_id, label) }
        .join(ch_mapping_grouped)
        .map { sample_id, label, haplotype_ids, results ->
            tuple(sample_id, label, haplotype_ids, results)
        }
        .set { ch_mapping_tagged }
    
    // Join everything with (sample_id, qc_label) as the key
    ch_quast_tagged
        .join(ch_merqury_tagged, by: [0, 1])  // Join on both sample_id AND qc_label
        .join(ch_busco_tagged, by: [0, 1])
        .join(ch_mapping_tagged, by: [0, 1])
        .map { sample_id, label, quast_results, merqury_results,
               haplotype_ids_busco, busco_results,
               haplotype_ids_mapping, mapping_results ->
            tuple(sample_id, label,
                  quast_results, merqury_results,
                  haplotype_ids_busco, busco_results,
                  haplotype_ids_mapping, mapping_results)
        }
        .set { ch_all_qc_final }
    
    COMBINE_ASSEMBLY_QC(ch_all_qc_final)
    
    emit:
    assembly_summary = COMBINE_ASSEMBLY_QC.out.summary
}