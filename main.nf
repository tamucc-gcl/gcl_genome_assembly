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
    - Optional decontamination at multiple stages
========================================================================================
*/

/*
========================================================================================
    PARAMETERS - FLATTENED FOR EASY COMMAND-LINE OVERRIDE
========================================================================================
    All parameters can now be overridden individually without affecting others.
    Example: --decon_run_on_contigs true --decon_source_taxid 373251
    
    IMPORTANT: Place these BEFORE any conditional logic that uses params
========================================================================================
*/

// ============================================================================
// Core pipeline parameters (keep existing ones)
// ============================================================================
params.sample_sheet = null
params.outdir = './results'
params.publish_dir_mode = 'link'

// Assembly QC parameters
params.busco_lineage = 'actinopterygii_odb10'
params.busco_downloads = '/work/birdlab/GCL/Databases/busco_datasets'

// Hi-C mapping and QC parameters
params.hic_coverage_window = 100000
params.hic_min_mapq = 30
params.hic_resolutions = "1000000,500000,100000,50000,10000"
params.hic_base_bin = "10000"
params.hic_plot_resolutions = "1000000,500000,100000"
params.hic_balance = false
params.hic_min_mapq_raw = 30
params.hic_min_mapq_filtered = 1

// Hi-C scaffolding parameters (YaHS)
params.yahs_min_contig_length = 10000
params.yahs_min_mapq = 1
params.yahs_resolutions = '10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000'
params.yahs_rounds_per_resolution = null
params.yahs_enzyme = null
params.yahs_no_contig_ec = false
params.yahs_no_scaffold_ec = false
params.bwa_mem2_hic_args = null
params.scaffold_min_size = 10000000

// Misassembly correction parameters
params.correct_contigs = true     // Run correction on contig assemblies
params.correct_scaffolds = true   // Run correction on scaffolded assemblies

// Contig correction parameters
params.contig_correction_min_depth = null  // null = use Inspector default (20% of average depth)
params.contig_correction_min_contig_length = 10000
params.contig_correction_min_contig_length_assemblyerror = 1000000
params.contig_correction_min_assembly_error_size = 50
params.contig_correction_max_assembly_error_size = 4000000

// Scaffold correction parameters (often need different thresholds)
params.scaffold_correction_min_depth = null  // null = use Inspector default (20% of average depth)
params.scaffold_correction_min_contig_length = 10000
params.scaffold_correction_min_contig_length_assemblyerror = 10000000  // Higher for scaffolds
params.scaffold_correction_min_assembly_error_size = 50
params.scaffold_correction_max_assembly_error_size = 10000000  // Higher for scaffolds

// ============================================================================
// Decontamination control (FLATTENED)
// ============================================================================
// When to run decontamination
params.decon_run_on_contigs = true       // Run on initial contig assemblies (before Hi-C mapping)
params.decon_run_on_scaffolds = false     // Run on scaffolded assemblies (after Hi-C scaffolding)

// Core settings
params.decon_source_taxid = 7898          // Actinopterygii; set your species taxid when possible

// Optional: adapter/vector removal
params.decon_run_fcs_adaptor = false      // Requires FCS-adaptor containers; enable only if configured
params.decon_fcsadaptor_mode = 'euk'      // 'euk' or 'prok'
params.decon_container_engine = 'singularity'

// Optional: generate evidence (coverage + taxonomy + blobtools plots)
params.decon_make_blobtools_evidence = true

// ============================================================================
// Database base directory
// ============================================================================
params.db_base = "/work/birdlab/databases"

// ============================================================================
// FCS-GX database settings (FLATTENED)
// ============================================================================
params.gxdb_dir = "${params.db_base}/fcs-gx"
params.gxdb_profile = 'all'          // 'all' | 'test-only'
params.gxdb_manifest = null          // optional: override manifest URL/path
params.gxdb_force = false            // re-download even if present

// ============================================================================
// DIAMOND / blobtools evidence DB settings (FLATTENED)
// ============================================================================
params.diamond_dmnd = null           // if you have a prebuilt .dmnd, set this
params.diamond_dir = "${params.db_base}/diamond"
params.diamond_name = "proteins"
params.diamond_taxdump_dir = "${params.db_base}/ncbi_taxonomy"
params.diamond_profile = 'custom'    // 'custom' (recommended)
params.diamond_fasta_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
params.diamond_taxonmap_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz"
params.diamond_force = false

// ============================================================================
// Evidence generation settings (FLATTENED)
// ============================================================================
params.evidence_map_preset = 'map-hifi'
params.evidence_diamond_max_target_seqs = 1
params.evidence_diamond_evalue = 1e-25
params.evidence_blob_min_contig_len = 1000

// ============================================================================
// BACKWARD COMPATIBILITY LAYER
// ============================================================================
// Reconstruct nested structures for internal workflow code
// This must come AFTER all the flattened parameters are defined above

params.gxdb = [
    dir: params.gxdb_dir,
    profile: params.gxdb_profile,
    manifest: params.gxdb_manifest,
    force: params.gxdb_force
]

params.diamond = [
    dmnd: params.diamond_dmnd,
    dir: params.diamond_dir,
    name: params.diamond_name,
    taxdump_dir: params.diamond_taxdump_dir,
    profile: params.diamond_profile,
    fasta_url: params.diamond_fasta_url,
    taxonmap_url: params.diamond_taxonmap_url,
    force: params.diamond_force
]

params.decon = [
    run_on_contigs: params.decon_run_on_contigs,
    run_on_scaffolds: params.decon_run_on_scaffolds,
    source_taxid: params.decon_source_taxid,
    run_fcs_adaptor: params.decon_run_fcs_adaptor,
    fcsadaptor_mode: params.decon_fcsadaptor_mode,
    container_engine: params.decon_container_engine,
    make_blobtools_evidence: params.decon_make_blobtools_evidence
]

params.evidence = [
    map_preset: params.evidence_map_preset,
    diamond_max_target_seqs: params.evidence_diamond_max_target_seqs,
    diamond_evalue: params.evidence_diamond_evalue,
    blob_min_contig_len: params.evidence_blob_min_contig_len
]

// Print pipeline header
log.info """\
    =========================================
    GENOME ASSEMBLY PIPELINE
    =========================================
    Sample sheet : ${params.sample_sheet}
    Output dir   : ${params.outdir}
    Misassembly correction: ${params.correct_contigs ? 'Contigs' : ''}${params.correct_scaffolds ? ' Scaffolds' : ''}${!params.correct_contigs && !params.correct_scaffolds ? 'Disabled' : ''}
    Decontamination: ${params.decon.run_on_contigs ? 'Contigs' : ''}${params.decon.run_on_scaffolds ? ' Scaffolds' : ''}${!params.decon.run_on_contigs && !params.decon.run_on_scaffolds ? 'Disabled' : ''}
    =========================================
    """
    .stripIndent()

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
include { ASSEMBLY_QC as ASSEMBLY_QC_CONTIG_CORRECTED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_CONTIG_DECONTAM } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_CORRECTED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_DECONTAM } from './workflows/assembly_qc.nf'

// HI-C MODULAR WORKFLOWS
include { HIC_QC_FROM_BAM as HIC_QC_FROM_BAM_RAW } from './workflows/hic_qc_from_bam.nf'
include { HIC_QC_FROM_BAM as HIC_QC_FROM_BAM_FILTERED } from './workflows/hic_qc_from_bam.nf'
include { HIC_SCAFFOLD_QC } from './workflows/hic_scaffold_qc.nf'

// DECONTAMINATION MODULAR WORKFLOWS
include { SETUP_DECONTAM_DBS } from './workflows/setup_decontam_dbs.nf'
include { DECONTAMINATE_ASSEMBLY as DECONTAMINATE_ASSEMBLY_CONTIG } from './workflows/decontaminate_assembly.nf'
include { DECONTAMINATE_ASSEMBLY as DECONTAMINATE_ASSEMBLY_SCAFFOLD } from './workflows/decontaminate_assembly.nf'
include { GENERATE_DECONTAM_EVIDENCE } from './workflows/generate_decontam_evidence.nf'

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
include { CORRECT_MISASSEMBLIES as CORRECT_MISASSEMBLIES_CONTIG } from './modules/correct_misassemblies.nf'
include { CORRECT_MISASSEMBLIES as CORRECT_MISASSEMBLIES_SCAFFOLD } from './modules/correct_misassemblies.nf'


/*
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
    ========================================================================================
        STEP 0: Setup Decontamination Databases (if enabled)
        Runs in parallel with BAM conversion and assembly
        Only executes if decontamination is requested
    ========================================================================================
    */
    if (params.decon.run_on_contigs || params.decon.run_on_scaffolds) {
        SETUP_DECONTAM_DBS()
        
        // Store outputs for later use
        ch_gxdb_dir = SETUP_DECONTAM_DBS.out.gxdb_dir
        ch_diamond_db = SETUP_DECONTAM_DBS.out.diamond_db
        ch_taxdump_dir = SETUP_DECONTAM_DBS.out.taxdump_dir
    }

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
        STEP 2: Trim Hi-C Reads
    ========================================================================================
    */
    TRIM_HIC(
        ch_input.map { sample_id, hifi_bam, hic_r1, hic_r2 ->
            tuple(sample_id, hic_r1, hic_r2)
        }
    )

    /*
    ========================================================================================
        STEP 3: Combine HiFi FASTQ with trimmed Hi-C reads
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
        STEP 4: Assemble with Hifiasm
    ========================================================================================
    */
    
    HIFIASM(ch_fastq_all)

    /*
    ========================================================================================
        STEP 4.5: Optional Misassembly Correction of Contig Assemblies
    ========================================================================================
    */
    if (params.correct_contigs) {
        // Prepare correction parameters for contig stage
        def contig_correction_params = [
            min_depth: params.contig_correction_min_depth,
            min_contig_length: params.contig_correction_min_contig_length,
            min_contig_length_assemblyerror: params.contig_correction_min_contig_length_assemblyerror,
            min_assembly_error_size: params.contig_correction_min_assembly_error_size,
            max_assembly_error_size: params.contig_correction_max_assembly_error_size
        ]

        // Split assemblies into individual haplotypes for correction
        HIFIASM.out.assemblies
            .flatMap { sample_id, hap1_fasta, hap2_fasta ->
                [
                    tuple("${sample_id}_hap1", hap1_fasta),
                    tuple("${sample_id}_hap2", hap2_fasta)
                ]
            }
            // make sample_id the FIRST element so it can join BAM_TO_FASTQ.out cleanly
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(sample_id, haplotype_id, fasta)
            }
            .join(BAM_TO_FASTQ.out)   // joins on sample_id (element 0)
            .map { sample_id, haplotype_id, fasta, hifi_reads ->
                tuple(haplotype_id, fasta, hifi_reads, "contig", contig_correction_params)
            }
            .set { ch_contigs_for_correction }
        
        ch_contigs_for_correction.view { "CONTIG_CORR_INPUT: $it" }
        CORRECT_MISASSEMBLIES_CONTIG(ch_contigs_for_correction)
        
        // Use corrected assemblies for downstream steps
        // Store them in channel format for both decontamination and Hi-C mapping
        CORRECT_MISASSEMBLIES_CONTIG.out.corrected
            .map { haplotype_id, fasta ->
                tuple(haplotype_id, fasta)
            }
            .set { ch_contigs_for_decontam_or_hic }
    } else {
        // Use original assemblies if correction not requested
        HIFIASM.out.assemblies
            .flatMap { sample_id, hap1_fasta, hap2_fasta ->
                [
                    tuple("${sample_id}_hap1", hap1_fasta),
                    tuple("${sample_id}_hap2", hap2_fasta)
                ]
            }
            .set { ch_contigs_for_decontam_or_hic }
    }

    /*
    ========================================================================================
        STEP 5: Optional Decontamination of Contig Assemblies
        Runs in parallel with Hi-C mapping preparation
        Databases were already set up in STEP 0 (parallel with assembly)
    ========================================================================================
    */
    if (params.decon.run_on_contigs) {
        // Run decontamination (parallel across all haplotypes)
        DECONTAMINATE_ASSEMBLY_CONTIG(
            ch_contigs_for_decontam_or_hic,
            ch_gxdb_dir
        )
        
        
        // Use decontaminated assemblies for downstream Hi-C mapping
        // Need to add sample_id back for joining with Hi-C reads
        DECONTAMINATE_ASSEMBLY_CONTIG.out.decontaminated
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(haplotype_id, sample_id, fasta)
            }
            .set { ch_individual_haplotypes }
    } else {
        // Use original assemblies if decontamination not requested
        ch_contigs_for_decontam_or_hic
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(haplotype_id, sample_id, fasta)
            }
            .set { ch_individual_haplotypes }
    }

    /*
    ========================================================================================
        STEP 6: Map Hi-C to Assemblies (contigs or decontaminated contigs)
    ========================================================================================
    */
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
        STEP 7: Filter Hi-C BAM Files
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
        STEP 8: Scaffold with Hi-C (PARALLEL WITH FILTERED CONTIG QC!)
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
        STEP 8.5: Optional Misassembly Correction of Scaffolded Assemblies
    ========================================================================================
    */
    if (params.correct_scaffolds) {
        // Prepare correction parameters for scaffold stage
        def scaffold_correction_params = [
            min_depth: params.scaffold_correction_min_depth,
            min_contig_length: params.scaffold_correction_min_contig_length,
            min_contig_length_assemblyerror: params.scaffold_correction_min_contig_length_assemblyerror,
            min_assembly_error_size: params.scaffold_correction_min_assembly_error_size,
            max_assembly_error_size: params.scaffold_correction_max_assembly_error_size
        ]

        // Correct scaffolds using HiFi reads
        SCAFFOLD_HIC.out.scaffolds
            .map { haplotype_id, scaffold ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(sample_id, haplotype_id, scaffold)
            }
            .join(BAM_TO_FASTQ.out)   // joins on sample_id (element 0)
            .map { sample_id, haplotype_id, scaffold, hifi_reads ->
                tuple(haplotype_id, scaffold, hifi_reads, "scaffold", scaffold_correction_params)
            }
            .set { ch_scaffolds_for_correction }
        
        ch_scaffolds_for_correction.view { "SCAFFOLD_CORR_INPUT: $it" }
        CORRECT_MISASSEMBLIES_SCAFFOLD(ch_scaffolds_for_correction)
        
        // Use corrected scaffolds for downstream steps
        CORRECT_MISASSEMBLIES_SCAFFOLD.out.corrected
            .map { haplotype_id, fasta ->
                tuple(haplotype_id, fasta)
            }
            .set { ch_scaffolds_for_decontam_or_final }
    } else {
        // Use original scaffolds if correction not requested
        SCAFFOLD_HIC.out.scaffolds
            .set { ch_scaffolds_for_decontam_or_final }
    }

    /*
    ========================================================================================
        STEP 9: Optional Decontamination of Scaffolded Assemblies
        Runs after scaffolding is complete
        Databases were already set up in STEP 0
    ========================================================================================
    */
    if (params.decon.run_on_scaffolds) {
        // Decontaminate scaffolds (parallel across all haplotypes)
        DECONTAMINATE_ASSEMBLY_SCAFFOLD(
            ch_scaffolds_for_decontam_or_final,
            ch_gxdb_dir
        )

        // Store final scaffolds
        ch_final_scaffolds = DECONTAMINATE_ASSEMBLY_SCAFFOLD.out.decontaminated
    } else {
        // Use scaffolds (corrected or original) without decontamination
        ch_final_scaffolds = ch_scaffolds_for_decontam_or_final
    }
    

    /*
    ========================================================================================
        STEP 16: Rescaffold broken assemblies (FUTURE)
    ========================================================================================
    */
    // 1. remap Hi-C reads to broken assemblies
    // 2. scaffold again with YaHS
    // 3. Assembly QC of rescaffolded assemblies
    /*
    ========================================================================================
        STEP 17: Gap Filling
    ========================================================================================
    */
    // 1. gap fill with HiFi reads using tgs gap closer
    // 2. Assembly QC of gap-filled assemblies
    // 3. contact-map based QC of gap-filled assemblies
    
    /*
    ========================================================================================
        STEP 18: Finalization
    ========================================================================================
    */
    // 1. Generate final QC reports
    //   - join all assembly QC summaries
    // 2. Dotplots of final assemblies vs each other (hap1 vs hap2)
    // 3. 

    /*
    ========================================================================================
        QC Steps
    ========================================================================================
    */

    /*
    ========================================================================================
        Sequencing QC
    ========================================================================================
    */
    /*
    ========================================================================================
        QC Raw Hi-C Reads
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
        QC HiFi Reads
    ========================================================================================
    */
    HIFI_QC(
        BAM_TO_FASTQ.out
    )

    /*
    ========================================================================================
        QC Trimmed Hi-C Reads
    ========================================================================================
    */
    HIC_QC_TRIMMED(
        TRIM_HIC.out.trimmed_reads,
        "trimmed"
    )
    
    /*
    ========================================================================================
        Assembly QC
    ========================================================================================
    */

    /*
    ========================================================================================
        QC Contig Genomes Assemblies
    ========================================================================================
    */
    ASSEMBLY_QC_INITIAL(
        HIFIASM.out.assemblies,
        BAM_TO_FASTQ.out,
        'contig'
    )

    /*
    ========================================================================================
        QC Contig misassembly corrected Genome Assemblies
    ========================================================================================
    */
    if (params.correct_contigs) {
        // QC on misassembly corrected assemblies
        // Re-pair haplotypes by sample for QC
        CORRECT_MISASSEMBLIES_CONTIG.out.corrected
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
                tuple(sample_id, hap_num, fasta)
            }
            .groupTuple()
            .map { sample_id, hap_nums, fastas ->
                // Sort by haplotype number to ensure hap1, hap2 order
                def sorted = [hap_nums, fastas].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_corrected_contig_paired }

        ASSEMBLY_QC_CONTIG_CORRECTED(
            ch_corrected_contig_paired,
            BAM_TO_FASTQ.out,
            'contig_corrected'
        )
    }

    /*
    ========================================================================================
        QC Contig Decontaminated Genome Assemblies
    ========================================================================================
    */
    if (params.decon.run_on_contigs) {
        // QC on decontaminated assemblies
        // Re-pair haplotypes by sample for QC
        DECONTAMINATE_ASSEMBLY_CONTIG.out.decontaminated
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
                tuple(sample_id, hap_num, fasta)
            }
            .groupTuple()
            .map { sample_id, hap_nums, fastas ->
                // Sort by haplotype number to ensure hap1, hap2 order
                def sorted = [hap_nums, fastas].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_decontam_paired }

        ASSEMBLY_QC_CONTIG_DECONTAM(
            ch_decontam_paired,
            BAM_TO_FASTQ.out,
            'contig_decontam'
        )
    }

    /*
    ========================================================================================
        QC Scaffolded Genomes
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
        QC Scaffold Misassembled Corrected Genome Assemblies
    ========================================================================================
    */
    if (params.correct_scaffolds) {
        // QC on misassembly corrected assemblies
        // Re-pair haplotypes by sample for QC
        CORRECT_MISASSEMBLIES_SCAFFOLD.out.corrected
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
                tuple(sample_id, hap_num, fasta)
            }
            .groupTuple()
            .map { sample_id, hap_nums, fastas ->
                // Sort by haplotype number to ensure hap1, hap2 order
                def sorted = [hap_nums, fastas].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_corrected_scaffold_paired }

        ASSEMBLY_QC_SCAFFOLD_CORRECTED(
            ch_corrected_scaffold_paired,
            BAM_TO_FASTQ.out,
            'scaffold_corrected'
        )
    }


    /*
    ========================================================================================
        QC Scaffold Decontaminated Genome Assemblies
    ========================================================================================
    */
    if (params.decon.run_on_scaffolds) {
        // QC on decontaminated assemblies
        // Re-pair haplotypes by sample for QC
        DECONTAMINATE_ASSEMBLY_SCAFFOLD.out.decontaminated
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
                tuple(sample_id, hap_num, fasta)
            }
            .groupTuple()
            .map { sample_id, hap_nums, fastas ->
                // Sort by haplotype number to ensure hap1, hap2 order
                def sorted = [hap_nums, fastas].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_decontam_scaffold_paired }

        ASSEMBLY_QC_SCAFFOLD_DECONTAM(
            ch_decontam_scaffold_paired,
            BAM_TO_FASTQ.out,
            'scaffold_decontam'
        )
    }

    /*
    ========================================================================================
        Hi-C Mapping QC
    ========================================================================================
    */
    /*
    ========================================================================================
        Hi-C Mapping QC on Raw BAMs
    ========================================================================================
    */
    //skip for development
    /*
    HIC_QC_FROM_BAM_RAW(
        MAP_HIC_TO_ASSEMBLY.out.bam,
        ch_assemblies_for_qc,
        "raw"
    )
    */
    /*
    ========================================================================================
         Hi-C Mapping QC on Filtered BAMs & contig (decontam if that option chosen) assembly
    ========================================================================================
    */
    //skip for development
    /*
    HIC_QC_FROM_BAM_FILTERED(
        FILTER_HIC_BAM.out.bam,
        ch_assemblies_for_qc,
        "filtered"
    )
    */
    /*
    ========================================================================================
        Hi-C Mapping QC on Scaffolds (MODULAR - LIFTOVER + QC)
    ========================================================================================
    */
    //skip for development
    /*
    // Use pairs.gz from filtered contig QC and lift to scaffold coordinates
    HIC_SCAFFOLD_QC(
        HIC_QC_FROM_BAM_FILTERED.out.pairs,
        SCAFFOLD_HIC.out.agp,
        SCAFFOLD_HIC.out.scaffolds,
        "filtered"
    )
    */
    /*
    ========================================================================================
        Generate Optional Decontamination Evidence
    ========================================================================================
    */
    if (params.decon.run_on_contigs) {
        // Optional: Generate evidence for decontamination decisions
        // This runs in parallel with Hi-C mapping preparation
        if (params.decon.make_blobtools_evidence) {
            GENERATE_DECONTAM_EVIDENCE(
                DECONTAMINATE_ASSEMBLY.out.decontaminated,
                DECONTAMINATE_ASSEMBLY.out.contaminants,
                DECONTAMINATE_ASSEMBLY.out.action_report,
                DECONTAMINATE_ASSEMBLY.out.taxonomy_report,
                BAM_TO_FASTQ.out,
                ch_diamond_db,
                ch_taxdump_dir
            )
        }
    }

  if (params.decon.run_on_scaffolds) {
        // Optional: Generate evidence for scaffold decontamination
        // This runs in parallel with scaffold QC
        if (params.decon.make_blobtools_evidence) {
            GENERATE_DECONTAM_EVIDENCE(
                DECONTAMINATE_ASSEMBLY.out.decontaminated,
                DECONTAMINATE_ASSEMBLY.out.contaminants,
                DECONTAMINATE_ASSEMBLY.out.action_report,
                DECONTAMINATE_ASSEMBLY.out.taxonomy_report,
                BAM_TO_FASTQ.out,
                ch_diamond_db,
                ch_taxdump_dir
            )
        }
  }
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