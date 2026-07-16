/*
========================================================================================
    DECONTAMINATE ASSEMBLY WORKFLOW
========================================================================================
    Repo location: workflows/decontaminate_assembly.nf

    - Screen assemblies for contaminants using NCBI FCS-GX
    - Optional: adapter/vector removal with FCS-adaptor
    - Clean genomes by removing identified contaminants
    - Per-haplotype throughout (each haplotype independent — critical for caching)
    - Reusable for contigs OR scaffolds via the `stage` string

    Channel topology is unchanged from the pre-meta version (the consume-twice pattern
    that runs in production); only the key type changed from haplotype_id to meta. The
    internal .join(by: 0) calls key on meta, which the FCS_* processes pass through unchanged.

    Usage:
    - DECONTAMINATE_ASSEMBLY(assemblies, gxdb_dir, "contig")
    - DECONTAMINATE_ASSEMBLY(assemblies, gxdb_dir, "scaffold")
========================================================================================
*/

nextflow.enable.dsl=2

include { FCS_ADAPTOR }      from '../modules/fcs_adaptor.nf'
include { FCS_GX_SCREEN }    from '../modules/fcs_gx_screen.nf'
include { FCS_CLEAN_GENOME } from '../modules/fcs_clean_genome.nf'

workflow DECONTAMINATE_ASSEMBLY {
    take:
    assemblies      // channel: tuple(meta, assembly_fasta)
    gxdb_dir        // channel: already a channel from SETUP_DECONTAM_DBS
    stage           // string: "contig" or "scaffold" for output directory organization

    main:

    /*
    ========================================================================================
        STEP 1: Optional Adapter/Vector Screening
    ========================================================================================
    */
    if (params.decon?.run_fcs_adaptor ?: false) {
        assemblies
            .map { meta, assembly_fasta ->
                tuple(meta, assembly_fasta, params.decon?.fcsadaptor_mode ?: 'euk', params.decon?.container_engine ?: 'singularity', stage)
            }
            .set { ch_adaptor_input }

        FCS_ADAPTOR(ch_adaptor_input)

        // Restore meta with cleaned assemblies (FCS_ADAPTOR passes meta through unchanged)
        assemblies
            .join(FCS_ADAPTOR.out.cleaned_fasta, by: 0)
            .map { meta, orig_fa, cleaned_fa -> tuple(meta, cleaned_fa) }
            .set { ch_cleaned_input }
    } else {
        // Skip adaptor step - use original assemblies
        assemblies.set { ch_cleaned_input }
    }

    /*
    ========================================================================================
        STEP 2: FCS-GX Screening (parallel across all haplotypes)
    ========================================================================================
    */
    // Per-sample source taxid: meta.taxid (sheet `taxid` column or params.taxid), falling
    // back to the global params.decon_source_taxid (nested: params.decon.source_taxid).
    ch_cleaned_input
        .map { meta, assembly_fasta ->
            def tax = meta.taxid ?: params.decon?.source_taxid
            if (tax == null)
                throw new IllegalArgumentException("Sample '${meta.id}': decontamination requires a taxid — set a per-row 'taxid' in the sample sheet, or --taxid / --decon_source_taxid.")
            tuple(meta, assembly_fasta, tax)
        }
        .combine(gxdb_dir)
        .combine(Channel.value(stage))
        .set { ch_fcs_gx_input }

    FCS_GX_SCREEN(ch_fcs_gx_input)

    /*
    ========================================================================================
        STEP 3: Clean Genome (remove contaminants)
        Join cleaned assemblies with their action reports by meta.
    ========================================================================================
    */
    ch_cleaned_input
        .join(FCS_GX_SCREEN.out.action_report, by: 0)
        .map { meta, assembly_fasta, action_report ->
            tuple(meta, assembly_fasta, action_report, stage)
        }
        .set { ch_clean_input }

    FCS_CLEAN_GENOME(ch_clean_input)

    emit:
    decontaminated = FCS_CLEAN_GENOME.out.decontaminated_fasta  // tuple(meta, fasta)
    contaminants = FCS_CLEAN_GENOME.out.contaminants_fasta      // tuple(meta, fasta)
    action_report = FCS_GX_SCREEN.out.action_report             // tuple(meta, report)
    taxonomy_report = FCS_GX_SCREEN.out.taxonomy_report         // tuple(meta, report)
}
