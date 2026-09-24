include { DECONTAMINATE_ASSEMBLY as DECONTAMINATE_ASSEMBLY_CONTIG } from '../workflows/decontaminate_assembly.nf'
include { FILTER_ORGANELLE } from '../modules/filter_organelle.nf'
include { REDUNDANS }      from '../modules/redundans.nf'
include { PILON }          from '../modules/pilon.nf'
include { PURGE_DUPS } from '../modules/purge_dups.nf'
include { CORRECT_MISASSEMBLIES as CORRECT_MISASSEMBLIES_CONTIG } from '../modules/correct_misassemblies.nf'

/* Contig conditioning, before mapping and scaffolding. */
workflow CONTIG_REFINEMENT {
    take:
    ch_contigs
    ch_organelle_baits
    ch_hifi_reads
    ch_shortread_reads
    ch_gxdb_dir

    capabilities

    main:
    ch_purged = Channel.empty()
    ch_shortread_conditioned = Channel.empty()
    ch_versions = Channel.empty()
    ch_corrected = Channel.empty()
    ch_decontaminated = Channel.empty()

    ch_contigs
        .map { meta, fasta -> tuple(meta.sample, meta, fasta) }
        .combine( ch_organelle_baits, by: 0 )
        .branch { sample, meta, fasta, baits ->
            has_bait: baits.size() > 0
            no_bait:  true
        }
        .set { ch_contigs_baitsplit }

    ch_contigs_baitsplit.has_bait
        .map { sample, meta, fasta, baits -> tuple(meta, fasta, baits) }
        .set { ch_organelle_filter_input }

    FILTER_ORGANELLE(ch_organelle_filter_input)

    ch_organelle_filtered = FILTER_ORGANELLE.out.filtered
        .mix( ch_contigs_baitsplit.no_bait.map { s, meta, fasta, b -> tuple(meta, fasta) } )

    // Fork the FILTERED contigs into the read-type-specific conditioning paths.
    ch_organelle_filtered
        .branch { meta, fasta ->
            hifi:      meta.assembler == 'hifiasm'
            shortread: true
        }
        .set { ch_filtered_by_type }

    ch_organelle_filtered_hifi = ch_filtered_by_type.hifi   // HiFi filtered -> PURGE_DUPS / QC

    /*
    ====================================================================================
        STEP 4b: PURGE DUPLICATES (optional)
        ch_organelle_filtered_hifi is per-haplotype; attach per-sample HiFi reads (combine
        by sample), one PURGE_DUPS over both haplotypes.
    ====================================================================================
    */
        ch_organelle_filtered_hifi
            .filter { meta, fasta -> meta.dedup == 'purge_dups' }
            .map { meta, fasta -> [ meta.sample, meta, fasta ] }
            .combine( ch_hifi_reads.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
            .map { sample, meta, fasta, hifi_reads -> tuple(meta, fasta, hifi_reads) }
            .set { ch_purge_dups_input }

        if (capabilities.purge) {
        PURGE_DUPS(ch_purge_dups_input)
        ch_purged = PURGE_DUPS.out.purged_assembly
        ch_versions = ch_versions.mix(PURGE_DUPS.out.versions)
        }
        ch_hifiasm_output = ch_purged
            .mix(ch_organelle_filtered_hifi.filter { meta, fasta -> meta.dedup != 'purge_dups' })

    // Short-read conditioning: REDUNDANS (reduce/scaffold/gap-close) -> optional Pilon.
    // Fed by the organelle-filtered short-read contigs (was ch_contigs_by_type.shortread).
    ch_filtered_by_type.shortread
        .map { meta, fasta -> [ meta.sample, meta, fasta ] }
        .combine( ch_shortread_reads.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
        .map { sample, meta, fasta, r1, r2 -> tuple(meta, fasta, r1, r2) }
        .set { ch_redundans_input }

    if (capabilities.spades) {
    REDUNDANS(ch_redundans_input)
    ch_versions = ch_versions.mix(REDUNDANS.out.versions)

    if (params.run_pilon) {
        REDUNDANS.out.assembly
            .map { meta, fasta -> [ meta.sample, meta, fasta ] }
            .combine( ch_shortread_reads.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
            .map { sample, meta, fasta, r1, r2 -> tuple(meta, fasta, r1, r2) }
            .set { ch_pilon_input }
        PILON(ch_pilon_input)
        ch_shortread_conditioned = PILON.out.assembly
        ch_versions = ch_versions.mix(PILON.out.versions)
    } else {
        ch_shortread_conditioned = REDUNDANS.out.assembly
    }
    }

    /*
    ====================================================================================
        STEP 4.5: Optional Misassembly Correction of Contig Assemblies (Inspector)
    ====================================================================================
    */
    if (capabilities.hifiasm && params.run_inspector_contigs) {
        ch_hifiasm_output
            .map { meta, fasta -> [ meta.sample, meta, fasta ] }
            .combine( ch_hifi_reads.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
            .map { sample, meta, fasta, hifi_fastq ->
                def correction_params = [
                    min_depth: params.inspector_contig_min_depth,
                    min_contig_length: params.inspector_contig_min_contig_bp,
                    min_contig_length_assemblyerror: params.inspector_contig_assemblyerror_min_bp,
                    min_assembly_error_size: params.inspector_contig_min_assembly_error_bp,
                    max_assembly_error_size: params.inspector_contig_max_assembly_error_bp,
                    skip_baseerror: !params.inspector_contig_base_error_check
                ]
                tuple(meta, fasta, hifi_fastq, "contig", correction_params)
            }
            .set { ch_correction_input }

        CORRECT_MISASSEMBLIES_CONTIG(ch_correction_input)
        ch_corrected = CORRECT_MISASSEMBLIES_CONTIG.out.corrected
        ch_hifi_conditioned = ch_corrected   // per-hap (meta, fasta)
        ch_versions = ch_versions.mix(CORRECT_MISASSEMBLIES_CONTIG.out.versions)
    } else {
        ch_hifi_conditioned = ch_hifiasm_output                           // per-hap (meta, fasta)
    }

    /*
    ====================================================================================
        STEP 5: Optional Decontamination of Contig Assemblies
        HiFi-conditioned (mito-filter / purge / correct) and short-read-conditioned
        (redundans / pilon) both flow through the SAME optional decontamination — FCS-GX is
        genome-based, so short-read assemblies are screened too (per-sample taxid via
        meta.taxid; see decontaminate_assembly.nf). Short-read gets NO Inspector correction
        (long-read-only). After decontam: HiFi continues into Hi-C scaffolding / the HiFi-only
        bypass; short-read is finished (nothing to scaffold with) → straight to FINALIZE.
    ====================================================================================
    */
    ch_assemblies_for_decontam = ch_hifi_conditioned.mix(ch_shortread_conditioned)

    if (params.run_decon_contigs) {
        DECONTAMINATE_ASSEMBLY_CONTIG(ch_assemblies_for_decontam, ch_gxdb_dir, "contig")
        ch_decontaminated = DECONTAMINATE_ASSEMBLY_CONTIG.out.decontaminated
        ch_decontaminated_contigs = ch_decontaminated
        ch_versions = ch_versions.mix(DECONTAMINATE_ASSEMBLY_CONTIG.out.versions)
    } else {
        ch_decontaminated_contigs = ch_assemblies_for_decontam
    }

    emit:
    assembly = ch_decontaminated_contigs
    organelle_filtered = ch_organelle_filtered
    purged = ch_purged
    shortread_conditioned = ch_shortread_conditioned
    corrected = ch_corrected
    decontaminated = ch_decontaminated
    versions = ch_versions
}
