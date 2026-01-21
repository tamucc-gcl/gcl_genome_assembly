#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    Genome Assembly and Scaffolding Pipeline
========================================================================================
    Author: Jason Selwyn
    Description: Pipeline for assembling and scaffolding genomes from HiFi and Hi-C reads
    
    MODULAR DESIGN:
    - Maximum parallelization
    - Reusable QC workflows as "functions"
    - Clear separation of concerns
========================================================================================
*/

// Print pipeline header
log.info """\
    =========================================
    GENOME ASSEMBLY PIPELINE
    =========================================
    Sample sheet : ${params.sample_sheet}
    Output dir   : ${params.outdir}
    =========================================
    """
    .stripIndent()

/*
========================================================================================
    PARAMETERS
========================================================================================
*/

params.sample_sheet = null
params.outdir = './results'
params.publish_dir_mode = 'link' //change to copy at end
params.busco_lineage = 'actinopterygii_odb10'
params.busco_downloads = '/work/birdlab/GCL/Databases/busco_datasets'
params.hic_coverage_window = 100000
params.hic_min_mapq = 30  // Minimum mapping quality for valid Hi-C pairs
params.hic_resolutions = "1000000,500000,100000,50000,10000"
params.hic_base_bin = "10000"
params.hic_plot_resolutions = "1000000,500000,100000"
params.hic_balance = false
params.hic_min_mapq_raw = 30
params.hic_min_mapq_filtered = 1
params.yahs_min_contig_length = 10000
params.yahs_min_mapq = 1
params.yahs_resolutions = '10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000'
params.yahs_rounds_per_resolution = null   // corresponds to -R if set
params.yahs_enzyme = null                  // corresponds to -e if set
params.yahs_no_contig_ec = true
params.yahs_no_scaffold_ec = true
params.bwa_mem2_hic_args = null

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!params.sample_sheet) {
    exit 1, "Sample sheet not specified! Please provide --sample_sheet"
}

/*
========================================================================================
    IMPORT FUNCTIONS
========================================================================================
*/
include { parseSampleSheet } from './functions/parse_sample_sheet.nf'

/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/
include { HIC_QC as HIC_QC_RAW } from './workflows/hic_qc.nf'
include { HIC_QC as HIC_QC_TRIMMED } from './workflows/hic_qc.nf'
include { HIFI_QC } from './workflows/hifi_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_INITIAL } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD } from './workflows/assembly_qc.nf'

// NEW MODULAR WORKFLOWS
include { HIC_QC_FROM_BAM as HIC_QC_FROM_BAM_RAW } from './workflows/hic_qc_from_bam.nf'
include { HIC_QC_FROM_BAM as HIC_QC_FROM_BAM_FILTERED } from './workflows/hic_qc_from_bam.nf'
include { HIC_SCAFFOLD_QC } from './workflows/hic_scaffold_qc.nf'

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { BAM_TO_FASTQ } from './modules/bam_to_fastq.nf'
include { TRIM_HIC } from './modules/trim_hic.nf'
include { HIFIASM } from './modules/hifiasm.nf'
include { MAP_HIC_TO_ASSEMBLY } from './modules/map_hic_to_assembly.nf'
include { FILTER_HIC_BAM } from './modules/filter_hic_bam.nf'
include { SCAFFOLD_HIC } from './modules/scaffold_hic.nf'


/*
include { QC_SCAFFOLDS } from './modules/qc_scaffolds.nf'
include { GAP_FILLING } from './modules/gap_filling.nf'
include { QC_FINAL } from './modules/qc_final.nf'
*/


/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Parse sample sheet and create input channel
    ch_input = parseSampleSheet(params.sample_sheet)

    /*
    // Debug: Print all channel contents
    ch_input.view { sample_id, hifi_bam, hic_r1, hic_r2 ->
        """
        ========================================
        Sample ID  : ${sample_id}
        HiFi BAM   : ${hifi_bam}
        Hi-C R1    : ${hic_r1}
        Hi-C R2    : ${hic_r2}
        ========================================
        """
    }
    */
    
    /*
    ========================================================================================
        STEP 1: Convert BAM to FASTQ
    ========================================================================================
    */
    
    BAM_TO_FASTQ(
        ch_input.map { sample_id, hifi_bam, hic_r1, hic_r2 ->
            tuple(sample_id, hifi_bam)
        }
    )

    /*
    ========================================================================================
        STEP 2: QC Raw Hi-C Reads
    ========================================================================================
    */
    HIC_QC_RAW(
        ch_input.map { sample_id, hifi_bam, hic_r1, hic_r2 ->
            tuple(sample_id, hic_r1, hic_r2)
        },
        "raw"
    )

    
    /*
    ========================================================================================
        STEP 3: Trim Hi-C Reads
    ========================================================================================
    */
    TRIM_HIC(
        ch_input.map { sample_id, hifi_bam, hic_r1, hic_r2 ->
            tuple(sample_id, hic_r1, hic_r2)
        }
    )
    

    /*
    ========================================================================================
        STEP 4: QC Trimmed Hi-C Reads
    ========================================================================================
    */
    HIC_QC_TRIMMED(
        TRIM_HIC.out.trimmed_reads,
        "trimmed"
    )
    
    /*
    ========================================================================================
        STEP 5: QC HiFi Reads - runs after converting BAM -> FASTQ
    ========================================================================================
    */
    HIFI_QC(
        BAM_TO_FASTQ.out
    )

    /*
    ========================================================================================
        STEP 6: Combine HiFi FASTQ with trimmed Hi-C reads
    ========================================================================================
    */
    TRIM_HIC.out.trimmed_reads
        .join(BAM_TO_FASTQ.out)
        .map { sample_id, hic_r1_trim, hic_r2_trim, hifi_fastq ->
            tuple(sample_id, hifi_fastq, hic_r1_trim, hic_r2_trim)
        }
        .set { ch_fastq_all }

    /*
    ========================================================================================
        STEP 7: Assemble with Hifiasm
    ========================================================================================
    */
    
    HIFIASM(ch_fastq_all)

    /*
    ========================================================================================
        STEP 8: QC Assembled Genomes
    ========================================================================================
    */
    ASSEMBLY_QC_INITIAL(
        HIFIASM.out.assemblies,
        BAM_TO_FASTQ.out,
        'contig'
    )

    /*
    ========================================================================================
        STEP 9: Map Hi-C to Contigs
    ========================================================================================
    */
    // Split assemblies into individual haplotypes
    HIFIASM.out.assemblies
        .flatMap { sample_id, hap1_fasta, hap2_fasta ->
            [
                tuple("${sample_id}_hap1", sample_id, hap1_fasta),
                tuple("${sample_id}_hap2", sample_id, hap2_fasta)
            ]
        }
        .set { ch_individual_haplotypes }
    
    // Combine each haplotype with its corresponding trimmed Hi-C reads
    ch_individual_haplotypes
        .map { haplotype_id, sample_id, fasta ->
            tuple(sample_id, haplotype_id, fasta)
        }
        .combine(TRIM_HIC.out.trimmed_reads, by: 0)
        .map { sample_id, haplotype_id, fasta, hic_r1, hic_r2 ->
            tuple(haplotype_id, fasta, hic_r1, hic_r2)
        }
        .set { ch_hic_mapping_input }
    
    // Map Hi-C reads to assemblies
    MAP_HIC_TO_ASSEMBLY(ch_hic_mapping_input)

    // Prepare assemblies channel for QC (haplotype_id, assembly_fasta)
    ch_individual_haplotypes
        .map { haplotype_id, sample_id, fasta ->
            tuple(haplotype_id, fasta)
        }
        .set { ch_assemblies_for_qc }

    /*
    ========================================================================================
        STEP 10: Hi-C Mapping QC on Raw BAMs (MODULAR - RUNS IN PARALLEL)
    ========================================================================================
    */
    
    HIC_QC_FROM_BAM_RAW(
        MAP_HIC_TO_ASSEMBLY.out.bam,
        ch_assemblies_for_qc,
        "raw"
    )
    
    /*
    ========================================================================================
        STEP 11: Filter Hi-C BAM Files (PARALLEL WITH RAW QC)
    ========================================================================================
    */
    
    // Combine BAM files with assemblies for filtering
    MAP_HIC_TO_ASSEMBLY.out.bam
        .join(ch_assemblies_for_qc)
        .set { ch_bam_with_assembly }
    
    // Filter BAM files to remove invalid pairs and duplicates
    FILTER_HIC_BAM(ch_bam_with_assembly)
    
    /*
    ========================================================================================
        STEP 12A: Hi-C Mapping QC on Filtered BAMs (MODULAR - PARALLEL WITH SCAFFOLDING!)
    ========================================================================================
    */
    
    HIC_QC_FROM_BAM_FILTERED(
        FILTER_HIC_BAM.out.bam,
        ch_assemblies_for_qc,
        "filtered"
    )
    
    /*
    ========================================================================================
        STEP 12B: Scaffold with Hi-C (PARALLEL WITH FILTERED CONTIG QC!)
    ========================================================================================
    */
    
    // Prepare input for scaffolding: (haplotype_id, filtered_bam, bai, assembly_fasta)
    FILTER_HIC_BAM.out.bam
        .join(ch_assemblies_for_qc)
        .set { ch_scaffolding_input }
    
    // Run Hi-C scaffolding
    SCAFFOLD_HIC(ch_scaffolding_input)

    /*
    ========================================================================================
        STEP 13: Hi-C Mapping QC on Scaffolds (MODULAR - LIFTOVER + QC)
    ========================================================================================
    */
    
    // Use pairs.gz from filtered contig QC and lift to scaffold coordinates
    HIC_SCAFFOLD_QC(
        HIC_QC_FROM_BAM_FILTERED.out.pairs,
        SCAFFOLD_HIC.out.agp,
        SCAFFOLD_HIC.out.scaffolds,
        "filtered"
    )

    /*
    ========================================================================================
        STEP 14: QC Scaffolded Genomes
    ========================================================================================
    */
    // Re-pair scaffolded haplotypes by sample
    SCAFFOLD_HIC.out.scaffolds
        .map { haplotype_id, scaffold ->
            // Extract sample_id and haplotype number
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
            tuple(sample_id, hap_num, scaffold)
        }
        .groupTuple()
        .map { sample_id, hap_nums, scaffolds ->
            // Sort by haplotype number to ensure hap1, hap2 order
            def sorted = [hap_nums, scaffolds].transpose().sort { it[0] }
            tuple(sample_id, sorted[0][1], sorted[1][1])
        }
        .set { ch_scaffolds_paired }

    ASSEMBLY_QC_SCAFFOLD(
        ch_scaffolds_paired,
        BAM_TO_FASTQ.out,
        'scaffold'
    )
    
    /*
    ========================================================================================
        FUTURE: Gap Filling and Finalization
    ========================================================================================
    */
    /*
    // Example of using modular QC on gap-filled assemblies:
    
    GAP_FILLING(ch_gapfill_input)
    
    // Reuse scaffold pairs for gap-filled assembly QC (no remapping!)
    HIC_QC_FROM_PAIRS(
        HIC_SCAFFOLD_QC.out.lifted_pairs,
        GAP_FILLING.out.gapfilled_assemblies,
        "gapfilled"
    )
    */
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """\
        Pipeline completed!
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Results   : ${params.outdir}
        """
        .stripIndent()
}