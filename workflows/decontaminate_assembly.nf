/*
========================================================================================
    DECONTAMINATE ASSEMBLY WORKFLOW
========================================================================================
    Purpose:
    - Screen assemblies for contaminants using NCBI FCS-GX
    - Optional: adapter/vector removal with FCS-adaptor
    - Clean genomes by removing identified contaminants
    - Fully parallelized across haplotypes
    
    Design:
    - Takes pre-prepared databases as input
    - Each haplotype processes independently
    - Returns both clean and contaminant sequences
    - Can be applied to contigs OR scaffolds
========================================================================================
*/

nextflow.enable.dsl=2

include { FCS_ADAPTOR }      from '../modules/fcs_adaptor.nf'
include { FCS_GX_SCREEN }    from '../modules/fcs_gx_screen.nf'
include { FCS_CLEAN_GENOME } from '../modules/fcs_clean_genome.nf'

workflow DECONTAMINATE_ASSEMBLY {
    take:
    assemblies      // channel: tuple(haplotype_id, assembly_fasta)
    gxdb_dir        // channel: already a channel from SETUP_DECONTAM_DBS
    
    main:
    
    /*
    ========================================================================================
        STEP 1: Optional Adapter/Vector Screening
    ========================================================================================
    */
    if (params.decon?.run_fcs_adaptor ?: false) {
        assemblies
            .map { haplotype_id, assembly_fasta ->
                tuple(assembly_fasta,
                      params.decon?.fcsadaptor_mode ?: 'euk',
                      params.decon?.container_engine ?: 'singularity')
            }
            .set { ch_adaptor_input }
        
        FCS_ADAPTOR(ch_adaptor_input)
        
        // Restore haplotype_id after adaptor cleaning
        assemblies
            .map { haplotype_id, assembly_fasta -> haplotype_id }
            .combine(FCS_ADAPTOR.out.cleaned_fasta)
            .set { ch_cleaned_input }
    } else {
        // Skip adaptor step - use original assemblies
        assemblies.set { ch_cleaned_input }
    }
    
    /*
    ========================================================================================
        STEP 2: FCS-GX Screening (Parallel Across All Haplotypes)
    ========================================================================================
    */
    // Extract just the assembly files for FCS_GX_SCREEN
    ch_cleaned_input
        .map { haplotype_id, assembly_fasta -> assembly_fasta }
        .set { ch_assembly_files }
    
    // Create source_taxid channel
    ch_source_taxid = Channel.value(params.decon?.source_taxid ?: 7898)
    
    // Call FCS_GX_SCREEN with separate inputs
    FCS_GX_SCREEN(
        ch_assembly_files,
        ch_source_taxid,
        gxdb_dir
    )
    
    /*
    ========================================================================================
        STEP 3: Clean Genome (Remove Identified Contaminants)
    ========================================================================================
    */
    FCS_CLEAN_GENOME(
        ch_assembly_files,
        FCS_GX_SCREEN.out.action_report
    )
    
    /*
    ========================================================================================
        STEP 4: Restore Haplotype IDs for All Outputs
    ========================================================================================
    */
    // Restore haplotype_id by position matching (since processes run in order)
    ch_cleaned_input
        .map { haplotype_id, assembly_fasta -> haplotype_id }
        .toList()
        .set { ch_haplotype_ids }
    
    FCS_CLEAN_GENOME.out.decontaminated_fasta
        .toList()
        .combine(ch_haplotype_ids)
        .flatMap { fastas, ids ->
            [fastas, ids].transpose().collect { fasta, id -> tuple(id, fasta) }
        }
        .set { ch_decontam_with_id }
    
    FCS_CLEAN_GENOME.out.contaminants_fasta
        .toList()
        .combine(ch_haplotype_ids)
        .flatMap { fastas, ids ->
            [fastas, ids].transpose().collect { fasta, id -> tuple(id, fasta) }
        }
        .set { ch_contam_with_id }
    
    FCS_GX_SCREEN.out.action_report
        .toList()
        .combine(ch_haplotype_ids)
        .flatMap { reports, ids ->
            [reports, ids].transpose().collect { report, id -> tuple(id, report) }
        }
        .set { ch_action_with_id }
    
    FCS_GX_SCREEN.out.taxonomy_report
        .toList()
        .combine(ch_haplotype_ids)
        .flatMap { reports, ids ->
            [reports, ids].transpose().collect { report, id -> tuple(id, report) }
        }
        .set { ch_taxonomy_with_id }
    
    emit:
    decontaminated = ch_decontam_with_id    // tuple(haplotype_id, clean_fasta)
    contaminants = ch_contam_with_id        // tuple(haplotype_id, contam_fasta)
    action_report = ch_action_with_id       // tuple(haplotype_id, action_report)
    taxonomy_report = ch_taxonomy_with_id   // tuple(haplotype_id, taxonomy_report)
    stdout_log = FCS_GX_SCREEN.out.stdout_log
}