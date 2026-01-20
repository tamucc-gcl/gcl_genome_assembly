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
    
    Can process either initial contigs or scaffolded assemblies with appropriate labeling
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
    qc_label     // string: "contigs" or "scaffolds"
    
    main:
    
    /*
    ========================================================================================
        Add QC label to assemblies for downstream processes
    ========================================================================================
    */
    assemblies
        .map { sample_id, hap1_fasta, hap2_fasta ->
            tuple(sample_id, hap1_fasta, hap2_fasta, qc_label)
        }
        .set { ch_assemblies_labeled }
    
    /*
    ========================================================================================
        QUAST - Run on both haplotypes per sample
    ========================================================================================
    */
    QUAST(ch_assemblies_labeled)
    
    /*
    ========================================================================================
        MERQURY - Run on both haplotypes per sample with HiFi reads
    ========================================================================================
    */
    MERQURY(
        ch_assemblies_labeled
            .map { sample_id, hap1_fasta, hap2_fasta, label ->
                tuple(sample_id, hap1_fasta, hap2_fasta, label)
            }
            .join(hifi_reads)
            .map { sample_id, hap1_fasta, hap2_fasta, label, hifi_fastq ->
                tuple(sample_id, hap1_fasta, hap2_fasta, hifi_fastq, label)
            }
    )
    
    /*
    ========================================================================================
        Split haplotypes for per-haplotype QC
    ========================================================================================
    */
    ch_assemblies_labeled
        .flatMap { sample_id, hap1_fasta, hap2_fasta, label ->
            [
                tuple("${sample_id}_hap1", sample_id, hap1_fasta, label),
                tuple("${sample_id}_hap2", sample_id, hap2_fasta, label)
            ]
        }
        .set { ch_individual_haplotypes }
    
    /*
    ========================================================================================
        BUSCO - Run on each haplotype independently
    ========================================================================================
    */
    BUSCO(
        ch_individual_haplotypes.map { haplotype_id, sample_id, fasta, label ->
            tuple(haplotype_id, fasta, label)
        }
    )
    
    /*
    ========================================================================================
        MAPPING QC - Map HiFi reads back to each haplotype
    ========================================================================================
    */
    // Combine each haplotype with its corresponding HiFi reads
    ch_individual_haplotypes
        .map { haplotype_id, sample_id, fasta, label ->
            tuple(sample_id, haplotype_id, fasta, label)
        }
        .combine(hifi_reads, by: 0)
        .map { sample_id, haplotype_id, fasta, label, hifi_fastq ->
            tuple(haplotype_id, fasta, hifi_fastq, label)
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
        .map { sample_id, results, label -> tuple(sample_id, label, results) }
        .set { ch_quast_by_sample }
    
    MERQURY.out.results
        .map { sample_id, results, label -> tuple(sample_id, label, results) }
        .set { ch_merqury_by_sample }
    
    BUSCO.out.results
        .map { haplotype_id, results, label -> 
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            tuple(sample_id, label, haplotype_id, results)
        }
        .groupTuple(by: [0, 1])
        .set { ch_busco_by_sample }
    
    MAPPING_QC.out.results
        .map { haplotype_id, results, label ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            tuple(sample_id, label, haplotype_id, results)
        }
        .groupTuple(by: [0, 1])
        .set { ch_mapping_by_sample }
    
    // Combine all QC results per sample and label
    ch_quast_by_sample
        .join(ch_merqury_by_sample, by: [0, 1])
        .join(ch_busco_by_sample, by: [0, 1])
        .join(ch_mapping_by_sample, by: [0, 1])
        .set { ch_all_qc }
    
    COMBINE_ASSEMBLY_QC(ch_all_qc)
    
    emit:
    assembly_summary = COMBINE_ASSEMBLY_QC.out.summary
}