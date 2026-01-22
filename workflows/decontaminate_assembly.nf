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
    ch_cleaned_input
        .combine(gxdb_dir)  // ✅ FIXED: No Channel.value() wrapper
        .map { haplotype_id, assembly_fasta, gxdb ->
            tuple(assembly_fasta,
                  params.decon?.source_taxid ?: 7898,
                  gxdb)
        }
        .set { ch_screen_input }
    
    FCS_GX_SCREEN(ch_screen_input)
    
    /*
    ========================================================================================
        STEP 3: Clean Genome (Remove Identified Contaminants)
    ========================================================================================
    */
    ch_cleaned_input
        .map { haplotype_id, assembly_fasta -> assembly_fasta }
        .combine(FCS_GX_SCREEN.out.action_report)
        .set { ch_clean_input }
    
    FCS_CLEAN_GENOME(ch_clean_input)
    
    /*
    ========================================================================================
        STEP 4: Restore Haplotype IDs for All Outputs
    ========================================================================================
    */
    // Decontaminated assemblies
    assemblies
        .map { haplotype_id, assembly_fasta -> haplotype_id }
        .combine(FCS_CLEAN_GENOME.out.decontaminated_fasta)
        .set { ch_decontam_with_id }
    
    // Contaminant sequences
    assemblies
        .map { haplotype_id, assembly_fasta -> haplotype_id }
        .combine(FCS_CLEAN_GENOME.out.contaminants_fasta)
        .set { ch_contam_with_id }
    
    // Action reports
    assemblies
        .map { haplotype_id, assembly_fasta -> haplotype_id }
        .combine(FCS_GX_SCREEN.out.action_report)
        .set { ch_action_with_id }
    
    // Taxonomy reports
    assemblies
        .map { haplotype_id, assembly_fasta -> haplotype_id }
        .combine(FCS_GX_SCREEN.out.taxonomy_report)
        .set { ch_taxonomy_with_id }
    
    emit:
    decontaminated = ch_decontam_with_id    // tuple(haplotype_id, clean_fasta)
    contaminants = ch_contam_with_id        // tuple(haplotype_id, contam_fasta)
    action_report = ch_action_with_id       // tuple(haplotype_id, action_report)
    taxonomy_report = ch_taxonomy_with_id   // tuple(haplotype_id, taxonomy_report)
    stdout_log = FCS_GX_SCREEN.out.stdout_log
}