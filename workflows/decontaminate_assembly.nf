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
        STEP 2: FCS-GX Screening (keep haplotype_id throughout)
    ========================================================================================
    */
    ch_cleaned_input
        .combine(Channel.value(params.decon?.source_taxid ?: 7898))
        .combine(gxdb_dir)
        .set { ch_gx_input }
    
    // We'll need a modified FCS_GX_SCREEN that accepts haplotype_id
    // For now, extract just files, then match back by position
    ch_cleaned_input
        .map { haplotype_id, assembly_fasta -> assembly_fasta }
        .set { ch_assembly_files }
    
    ch_source_taxid = Channel.value(params.decon?.source_taxid ?: 7898)
    
    FCS_GX_SCREEN(
        ch_assembly_files,
        ch_source_taxid,
        gxdb_dir
    )
    
    /*
    ========================================================================================
        STEP 3: Restore haplotype IDs and prepare for cleaning
    ========================================================================================
    */
    ch_cleaned_input
        .map { haplotype_id, assembly_fasta -> haplotype_id }
        .collect()
        .set { ch_haplotype_ids_list }
    
    // Match action reports with assemblies by position
    ch_cleaned_input
        .collect()
        .combine(FCS_GX_SCREEN.out.action_report.collect())
        .flatMap { assemblies, reports ->
            assemblies.withIndex().collect { asm, idx ->
                tuple(asm[0], asm[1], reports[idx])  // (haplotype_id, assembly, action_report)
            }
        }
        .set { ch_clean_input }
    
    FCS_CLEAN_GENOME(ch_clean_input)
    
    // Match taxonomy reports
    FCS_GX_SCREEN.out.taxonomy_report
        .collect()
        .combine(ch_haplotype_ids_list)
        .flatMap { reports, ids ->
            reports.withIndex().collect { report, idx -> tuple(ids[idx], report) }
        }
        .set { ch_taxonomy_with_id }
    
    // Match action reports
    FCS_GX_SCREEN.out.action_report
        .collect()
        .combine(ch_haplotype_ids_list)
        .flatMap { reports, ids ->
            reports.withIndex().collect { report, idx -> tuple(ids[idx], report) }
        }
        .set { ch_action_with_id }
    
    emit:
    decontaminated = FCS_CLEAN_GENOME.out.decontaminated_fasta  // tuple(haplotype_id, fasta)
    contaminants = FCS_CLEAN_GENOME.out.contaminants_fasta      // tuple(haplotype_id, fasta)
    action_report = ch_action_with_id                           // tuple(haplotype_id, report)
    taxonomy_report = ch_taxonomy_with_id                       // tuple(haplotype_id, report)
    stdout_log = FCS_GX_SCREEN.out.stdout_log
}