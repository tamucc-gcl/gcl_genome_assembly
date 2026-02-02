/*
========================================================================================
    DECONTAMINATE ASSEMBLY WORKFLOW (FIXED FOR CACHING)
========================================================================================
    Purpose:
    - Screen assemblies for contaminants using NCBI FCS-GX
    - Optional: adapter/vector removal with FCS-adaptor
    - Clean genomes by removing identified contaminants
    - Fully parallelized across haplotypes with proper sample independence
    
    Design:
    - Takes pre-prepared databases as input
    - Each haplotype processes independently (CRITICAL FOR CACHING)
    - Returns both clean and contaminant sequences
    - Can be applied to contigs OR scaffolds (controlled by stage parameter)
    - Stage parameter controls output directory organization
    
    CACHING FIX:
    - Maintains haplotype_id throughout the workflow
    - Each sample processed completely independently
    - Adding new samples won't invalidate cache for existing samples
    
    Usage:
    - DECONTAMINATE_ASSEMBLY(assemblies, gxdb_dir, "contig")  // For contig decontamination
    - DECONTAMINATE_ASSEMBLY(assemblies, gxdb_dir, "scaffold") // For scaffold decontamination
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
    stage           // string: "contig" or "scaffold" for output directory organization
    
    main:
    
    /*
    ========================================================================================
        STEP 1: Optional Adapter/Vector Screening
        CACHING FIX: Process each haplotype independently
    ========================================================================================
    */
    if (params.decon?.run_fcs_adaptor ?: false) {
        // Process each assembly file independently - maintain haplotype_id
        assemblies
            .map { haplotype_id, assembly_fasta ->
                tuple(haplotype_id, assembly_fasta, params.decon?.fcsadaptor_mode ?: 'euk', params.decon?.container_engine ?: 'singularity', stage)
            }
            .set { ch_adaptor_input }
        
        FCS_ADAPTOR(ch_adaptor_input)
        
        // Restore haplotype_id with cleaned assemblies
        assemblies
            .map { haplotype_id, assembly_fasta -> tuple(haplotype_id, assembly_fasta) }
            .join(
                FCS_ADAPTOR.out.cleaned_fasta.map { haplotype_id, cleaned -> tuple(haplotype_id, cleaned) }
            )
            .map { haplotype_id, orig_fa, cleaned_fa -> tuple(haplotype_id, cleaned_fa) }
            .set { ch_cleaned_input }
    } else {
        // Skip adaptor step - use original assemblies
        assemblies.set { ch_cleaned_input }
    }
    
    /*
    ========================================================================================
        STEP 2: FCS-GX Screening (Parallel Across All Haplotypes)
        CACHING FIX: Each haplotype processes completely independently
    ========================================================================================
    */
    // Prepare input for FCS_GX_SCREEN - each haplotype gets its own process invocation
    ch_cleaned_input
        .combine(Channel.value(params.decon?.source_taxid ?: 7898))
        .combine(gxdb_dir)
        .combine(Channel.value(stage))
        .set { ch_fcs_gx_input }
    
    // Each assembly screened independently - critical for caching
    FCS_GX_SCREEN(ch_fcs_gx_input)
    
    /*
    ========================================================================================
        STEP 3: Clean Genome (Remove Contaminants)
        CACHING FIX: Join by haplotype_id to maintain independence
    ========================================================================================
    */
    // Join cleaned assemblies with their corresponding action reports by haplotype_id
    ch_cleaned_input
        .join(FCS_GX_SCREEN.out.action_report, by: 0)
        .map { haplotype_id, assembly_fasta, action_report ->
            tuple(haplotype_id, assembly_fasta, action_report, stage)
        }
        .set { ch_clean_input }
    
    FCS_CLEAN_GENOME(ch_clean_input)
    
    emit:
    decontaminated = FCS_CLEAN_GENOME.out.decontaminated_fasta  // tuple(haplotype_id, fasta)
    contaminants = FCS_CLEAN_GENOME.out.contaminants_fasta      // tuple(haplotype_id, fasta)
    action_report = FCS_GX_SCREEN.out.action_report             // tuple(haplotype_id, report)
    taxonomy_report = FCS_GX_SCREEN.out.taxonomy_report         // tuple(haplotype_id, report)
}