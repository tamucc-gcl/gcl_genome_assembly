/*
========================================================================================
    ASSEMBLY QC SUBWORKFLOW
========================================================================================
    Comprehensive QC for phased genome assemblies
    - QUAST: per-sample (both haplotypes)
    - MERQURY: per-sample (both haplotypes)
    - BUSCO: per-haplotype
    - MAPPING: per-haplotype (map HiFi reads back to assembly)
    - COMBINE_QC: aggregate and visualize all QC metrics
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
    qc_label    // value/channel: label for output subfolder (e.g. 'contig' or 'scaffold')
    
    main:
    
    /*
    ========================================================================================
        QUAST - Run on both haplotypes per sample
    ========================================================================================
    */
    QUAST(assemblies)
    
    /*
    ========================================================================================
        MERQURY - Run on both haplotypes per sample with HiFi reads
    ========================================================================================
    */
    MERQURY(
        assemblies.join(hifi_reads)
    )
    
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
    // Combine each haplotype with its corresponding HiFi reads
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
    ========================================================================================
    */
    // Collect all QC outputs by sample
    QUAST.out.results
        .map { sample_id, results -> tuple(sample_id, results) }
        .set { ch_quast_by_sample }
    
    MERQURY.out.results
        .map { sample_id, results -> tuple(sample_id, results) }
        .set { ch_merqury_by_sample }
    
    BUSCO.out.results
        .map { haplotype_id, results -> 
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            tuple(sample_id, haplotype_id, results)
        }
        .groupTuple()
        .set { ch_busco_by_sample }
    
    MAPPING_QC.out.results
        .map { haplotype_id, results ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            tuple(sample_id, haplotype_id, results)
        }
        .groupTuple()
        .set { ch_mapping_by_sample }
    
    // Combine all QC results per sample
    ch_quast_by_sample
        .join(ch_merqury_by_sample)
        .join(ch_busco_by_sample)
        .join(ch_mapping_by_sample)
        .set { ch_all_qc }
    
    // Attach qc_label to each sample so outputs go to separate folders
    assemblies
        .map { sample_id, hap1_fasta, hap2_fasta -> tuple(sample_id, qc_label) }
        .set { ch_qc_label_by_sample }

    ch_all_qc
        .join(ch_qc_label_by_sample)
        .map { sample_id, quast_results, merqury_results, haplotype_ids_busco, busco_results, haplotype_ids_mapping, mapping_results, qc_label ->
            tuple(sample_id, qc_label,
                  quast_results,
                  merqury_results,
                  haplotype_ids_busco,
                  busco_results,
                  haplotype_ids_mapping,
                  mapping_results)
        }
        .set { ch_all_qc_labeled }
    COMBINE_ASSEMBLY_QC(ch_all_qc_labeled)
    
    emit:
    //quast_results = QUAST.out.results
    //merqury_results = MERQURY.out.results
    //busco_results = BUSCO.out.results
    //mapping_results = MAPPING_QC.out.results
    //combined_report = COMBINE_ASSEMBLY_QC.out.report
    //combined_plots = COMBINE_ASSEMBLY_QC.out.plots
    assembly_summary = COMBINE_ASSEMBLY_QC.out.summary
}