#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    Genome Assembly and Scaffolding Pipeline
========================================================================================
    Author: Your Name
    Description: Pipeline for assembling and scaffolding genomes from HiFi and Hi-C reads
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
params.publish_dir_mode = 'copy'

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
    IMPORT MODULES
========================================================================================
*/

include { parseSampleSheet } from './functions/parse_sample_sheet.nf'

/*
include { BAM_TO_FASTQ } from './modules/bam_to_fastq.nf'
include { FASTQC_HIC } from './modules/fastqc_hic.nf'
include { MULTIQC_HIC } from './modules/multiqc_hic.nf'
include { FASTQC_HIFI } from './modules/fastqc_hifi.nf'
include { MULTIQC_HIFI } from './modules/multiqc_hifi.nf'

include { HIFIASM } from './modules/hifiasm.nf'
include { QC_ASSEMBLY } from './modules/qc_assembly.nf'
include { SCAFFOLD_HIC } from './modules/scaffold_hic.nf'
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
    
    /*
    ========================================================================================
        STEP 1: Convert BAM to FASTQ
    ========================================================================================
    */
    /*
    BAM_TO_FASTQ(
        ch_input.map { sample_id, hifi_bam, hic_r1, hic_r2 ->
            tuple(sample_id, hifi_bam)
        }
    )
    
    // Combine HiFi FASTQ with Hi-C reads
    ch_input
        .map { sample_id, hifi_bam, hic_r1, hic_r2 ->
            tuple(sample_id, hic_r1, hic_r2)
        }
        .join(BAM_TO_FASTQ.out)
        .map { sample_id, hic_r1, hic_r2, hifi_fastq ->
            tuple(sample_id, hifi_fastq, hic_r1, hic_r2)
        }
        .set { ch_fastq_all }
    
    /*

    /*
    ========================================================================================
        STEP 2: QC Raw Hi-C Reads
    ========================================================================================
    */
    /*
    FASTQC_HIC(
        ch_fastq_all.map { sample_id, hifi_fastq, hic_r1, hic_r2 ->
            tuple(sample_id, hic_r1, hic_r2)
        }
    )
    
    MULTIQC_HIC(
        FASTQC_HIC.out.collect()
    )
    */
    /*
    ========================================================================================
        STEP 3: QC HiFi Reads
    ========================================================================================
    */
    /*
    FASTQC_HIFI(
        ch_fastq_all.map { sample_id, hifi_fastq, hic_r1, hic_r2 ->
            tuple(sample_id, hifi_fastq)
        }
    )
    
    MULTIQC_HIFI(
        FASTQC_HIFI.out.collect()
    )
    */
    /*
    ========================================================================================
        STEP 4: Assemble with Hifiasm
    ========================================================================================
    */
    /*
    HIFIASM(ch_fastq_all)
    */
    /*
    ========================================================================================
        STEP 5: QC Assembled Genomes
    ========================================================================================
    */
    /*
    // Hifiasm outputs phased assemblies (hap1 and hap2)
    // Flatten to run QC on each haplotype independently
    HIFIASM.out
        .flatMap { sample_id, hap1_fasta, hap2_fasta ->
            [
                tuple("${sample_id}_hap1", hap1_fasta),
                tuple("${sample_id}_hap2", hap2_fasta)
            ]
        }
        .set { ch_assemblies }
    
    QC_ASSEMBLY(ch_assemblies)
    */
    /*
    ========================================================================================
        STEP 6: Scaffold with Hi-C
    ========================================================================================
    */
    /*
    // Combine assemblies with Hi-C reads for scaffolding
    HIFIASM.out
        .join(
            ch_fastq_all.map { sample_id, hifi_fastq, hic_r1, hic_r2 ->
                tuple(sample_id, hic_r1, hic_r2)
            }
        )
        .flatMap { sample_id, hap1_fasta, hap2_fasta, hic_r1, hic_r2 ->
            [
                tuple("${sample_id}_hap1", hap1_fasta, hic_r1, hic_r2),
                tuple("${sample_id}_hap2", hap2_fasta, hic_r1, hic_r2)
            ]
        }
        .set { ch_scaffold_input }
    
    SCAFFOLD_HIC(ch_scaffold_input)
    */
    /*
    ========================================================================================
        STEP 7: QC Scaffolded Genomes
    ========================================================================================
    */
    /*
    QC_SCAFFOLDS(SCAFFOLD_HIC.out)
    */
    /*
    ========================================================================================
        STEP 8: Gap Filling and Finalization
    ========================================================================================
    */
    /*
    // Combine scaffolds with HiFi reads for gap filling
    SCAFFOLD_HIC.out
        .map { haplotype_id, scaffold ->
            // Extract original sample_id from haplotype_id (remove _hap1/_hap2)
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            tuple(sample_id, haplotype_id, scaffold)
        }
        .combine(
            ch_fastq_all.map { sample_id, hifi_fastq, hic_r1, hic_r2 ->
                tuple(sample_id, hifi_fastq)
            },
            by: 0
        )
        .map { sample_id, haplotype_id, scaffold, hifi_fastq ->
            tuple(haplotype_id, scaffold, hifi_fastq)
        }
        .set { ch_gapfill_input }
    
    GAP_FILLING(ch_gapfill_input)
    */
    /*
    ========================================================================================
        STEP 9: Final QC
    ========================================================================================
    */
    /*
    QC_FINAL(GAP_FILLING.out)
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