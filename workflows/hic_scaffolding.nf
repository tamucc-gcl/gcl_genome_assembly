include { DECONTAMINATE_ASSEMBLY as DECONTAMINATE_ASSEMBLY_SCAFFOLD } from '../workflows/decontaminate_assembly.nf'
include { CORRECT_MISASSEMBLIES as CORRECT_MISASSEMBLIES_SCAFFOLD } from '../modules/correct_misassemblies.nf'
include { MAP_HIC_TO_ASSEMBLY } from '../modules/map_hic_to_assembly.nf'
include { MAP_HIC_TO_ASSEMBLY as MAP_HIC_TO_SCAFFOLD } from '../modules/map_hic_to_assembly.nf'
include { FILTER_HIC_BAM } from '../modules/filter_hic_bam.nf'
include { FILTER_HIC_BAM as FILTER_HIC_BAM_SCAFFOLD } from '../modules/filter_hic_bam.nf'
include { SCAFFOLD_HIC as SCAFFOLD_HIC_ROUND1 } from '../modules/scaffold_hic.nf'
include { SCAFFOLD_HIC as SCAFFOLD_HIC_ROUND2 } from '../modules/scaffold_hic.nf'
include { GAP_FILLING } from '../modules/gap_filling.nf'
include { TELOCLIP_EXTEND; COLLECT_TELOCLIP_STATS } from '../modules/teloclip.nf'
include { HIC_BAM_METRICS as HIC_BAM_METRICS_CONTIG; HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_CONTIG } from '../modules/hic_mapping_metrics.nf'
include { HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_CONTIGSCAF } from '../modules/hic_mapping_metrics.nf'
include { HIC_BAM_METRICS as HIC_BAM_METRICS_SCAFFOLD; HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_SCAFFOLD } from '../modules/hic_mapping_metrics.nf'
include { HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_SCAFFOLDSCAF } from '../modules/hic_mapping_metrics.nf'

/* Hi-C mapping, scaffolding rounds and long-read finishing.
 * Sequence-editing operations and their QC outputs stay in this workflow.
 * Existing task aliases are retained within this permanent boundary.
 */
workflow HIC_SCAFFOLDING {
    take:
    ch_decontaminated_contigs
    ch_hifi_reads
    ch_hic_reads
    ch_telo_by_taxid
    ch_gxdb_dir

    capabilities

    main:
    ch_round1 = Channel.empty()
    ch_round1_agp = Channel.empty()
    ch_contig_pairs = Channel.empty()
    ch_filled = Channel.empty()
    ch_final_scaffolds_round2 = Channel.empty()
    ch_scaffold_round2_agp = Channel.empty()
    ch_all_bam_metrics = Channel.empty()
    ch_all_pairs_metrics = Channel.empty()
    ch_versions = Channel.empty()
    ch_scaffold_corrected = Channel.empty()
    ch_scaffold_decontaminated = Channel.empty()
    ch_extended = Channel.empty()

    ch_individual_haplotypes = ch_decontaminated_contigs.filter { meta, fasta -> meta.assembler == 'hifiasm' }
    ch_shortread_finished    = ch_decontaminated_contigs.filter { meta, fasta -> meta.assembler == 'spades' }

    // HiFi-only rows have no Hi-C, so they drop out of the Hi-C scaffolding block below.
    // Their decontaminated contigs are their "scaffolds" — rejoin at gap-filling.
    ch_hifi_only_scaffolds = ch_individual_haplotypes.filter { meta, fasta -> !meta.hic }

    if (capabilities.scaffold) {
    /*
    ====================================================================================
        STEP 6: Map Hi-C to Assemblies (contigs or decontaminated contigs)
    ====================================================================================
    */
    // Combine each haplotype with its sample's trimmed Hi-C reads (key on meta.sample —
    // TRIM_HIC carries the sample-level meta; ch_individual_haplotypes is per-haplotype)
    ch_individual_haplotypes
        .map { meta, fasta -> [ meta.sample, meta, fasta ] }
        .combine( ch_hic_reads.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
        .map { sample, meta, fasta, hic_r1, hic_r2 ->
            tuple(meta, fasta, hic_r1, hic_r2, "contig")
        }
        .set { ch_hic_mapping_input }

    MAP_HIC_TO_ASSEMBLY(ch_hic_mapping_input)
    ch_versions = ch_versions.mix(MAP_HIC_TO_ASSEMBLY.out.versions)

    // checkpoint: contig_raw_map (BAM-level only)
    MAP_HIC_TO_ASSEMBLY.out.bam
        .map { meta, stage, bam, bai -> tuple(meta, "contig_raw_map", bam, bai) }
        .set { ch_hic_raw_bam_for_qc }

    // Hi-C QC-metric accumulators — appended at each metrics step below and
    // consumed by COMPILE_FINAL_QC (kept here so the QC phase needs only these two channels).
    ch_all_bam_metrics   = Channel.empty()
    ch_all_pairs_metrics = Channel.empty()

    HIC_BAM_METRICS_CONTIG(ch_hic_raw_bam_for_qc)
    ch_all_bam_metrics = ch_all_bam_metrics.mix(HIC_BAM_METRICS_CONTIG.out.metrics)

    // Assemblies channel for the filter / scaffold joins (meta, fasta)
    ch_individual_haplotypes
        .map { meta, fasta -> tuple(meta, fasta) }
        .set { ch_assemblies_for_qc }


    /*
    ====================================================================================
        STEP 7: Filter Hi-C BAM Files
    ====================================================================================
    */
    MAP_HIC_TO_ASSEMBLY.out.bam
        .join(ch_assemblies_for_qc)
        .map { meta, stage, bam, bai, assembly_fasta ->
            tuple(meta, stage, bam, bai, assembly_fasta)
        }
        .set { ch_bam_with_assembly }

    FILTER_HIC_BAM(ch_bam_with_assembly)
    ch_versions = ch_versions.mix(FILTER_HIC_BAM.out.versions)

    // checkpoint: contig_filtered (pairs-level + retention); no AGP here -> []
    FILTER_HIC_BAM.out.pairs
        .join(FILTER_HIC_BAM.out.parse_stats)
        .join(FILTER_HIC_BAM.out.dedup_stats)
        .map { meta, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats ->
            tuple(meta, "contig_filtered", pairs_gz, [], parse_stats, dedup_stats)
        }
        .set { ch_hic_pairs_contig_filtered_for_qc }

    HIC_PAIRS_METRICS_CONTIG(ch_hic_pairs_contig_filtered_for_qc)
    ch_all_pairs_metrics = ch_all_pairs_metrics.mix(HIC_PAIRS_METRICS_CONTIG.out.metrics)
    
    /*
    ====================================================================================
        STEP 8: First Round of Scaffolding with Hi-C
    ====================================================================================
    */
    FILTER_HIC_BAM.out.bam
        .join(ch_assemblies_for_qc)
        .map { meta, stage, bam, bai, assembly_fasta ->
            tuple(meta, bam, bai, assembly_fasta, "round1", [ min_contig_length: params.yahs_round1_min_contig_bp, min_mapq: params.yahs_round1_min_mapq, resolutions: params.yahs_round1_resolutions, rounds_per_resolution: params.yahs_round1_rounds_per_resolution, enzyme: params.yahs_round1_enzyme, no_contig_ec: !params.yahs_round1_contig_ec, no_scaffold_ec: !params.yahs_round1_scaffold_ec ])
        }
        .set { ch_scaffolding_round1_input }

    SCAFFOLD_HIC_ROUND1(ch_scaffolding_round1_input)
    ch_versions = ch_versions.mix(SCAFFOLD_HIC_ROUND1.out.versions)

    // checkpoint: scaffold_space (relabel contigs->scaffolds via round1 AGP; no remap)
    FILTER_HIC_BAM.out.pairs
        .join(SCAFFOLD_HIC_ROUND1.out.agp)
        .join(FILTER_HIC_BAM.out.parse_stats)
        .join(FILTER_HIC_BAM.out.dedup_stats)
        .map { meta, stage, pairs_gz, agp, stage2, parse_stats, stage3, dedup_stats ->
            tuple(meta, "scaffold_space", pairs_gz, agp, parse_stats, dedup_stats)
        }
        .set { ch_hic_pairs_scaffold_space_for_qc }

    HIC_PAIRS_METRICS_CONTIGSCAF(ch_hic_pairs_scaffold_space_for_qc)
    ch_all_pairs_metrics = ch_all_pairs_metrics.mix(HIC_PAIRS_METRICS_CONTIGSCAF.out.metrics)

    /*
    ====================================================================================
        STEP 8.5: Optional Misassembly Correction of Scaffolded Assemblies (Inspector)
    ====================================================================================
    */
    if (params.run_inspector_scaffolds) {
        // Combine each scaffolded haplotype with its sample's HiFi reads (key on meta.sample)
        SCAFFOLD_HIC_ROUND1.out.scaffolds
            .map { meta, scaffold -> [ meta.sample, meta, scaffold ] }
            .combine( ch_hifi_reads.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
            .map { sample, meta, scaffold, hifi_fastq ->
                def correction_params = [
                    min_depth: params.inspector_scaffold_min_depth,
                    min_contig_length: params.inspector_scaffold_min_contig_bp,
                    min_contig_length_assemblyerror: params.inspector_scaffold_assemblyerror_min_bp,
                    min_assembly_error_size: params.inspector_scaffold_min_assembly_error_bp,
                    max_assembly_error_size: params.inspector_scaffold_max_assembly_error_bp,
                    skip_baseerror: !params.inspector_scaffold_base_error_check
                ]
                tuple(meta, scaffold, hifi_fastq, "scaffold", correction_params)
            }
            .set { ch_scaffold_correction_input }

        CORRECT_MISASSEMBLIES_SCAFFOLD(ch_scaffold_correction_input)
        ch_scaffold_corrected = CORRECT_MISASSEMBLIES_SCAFFOLD.out.corrected
        ch_scaffolds_for_decontam = ch_scaffold_corrected   // per-hap (meta, fasta)
    } else {
        ch_scaffolds_for_decontam = SCAFFOLD_HIC_ROUND1.out.scaffolds              // per-hap (meta, fasta)
    }

    /*
    ========================================================================================
        STEP 9: Optional Decontamination of Scaffolded Assemblies
        Runs after scaffolding (and optional correction) is complete
        Databases were already set up in STEP 0
        Works on either original scaffolds OR corrected scaffolds (if Inspector was run)
    ========================================================================================
    */
    if (params.run_decon_scaffolds) {
        // Decontaminate scaffolds (parallel across all haplotypes)
        DECONTAMINATE_ASSEMBLY_SCAFFOLD(
            ch_scaffolds_for_decontam,
            ch_gxdb_dir,
            "scaffold"
        )

        // Store final decontaminated scaffolds
        ch_scaffold_decontaminated = DECONTAMINATE_ASSEMBLY_SCAFFOLD.out.decontaminated
        ch_final_scaffolds = ch_scaffold_decontaminated
    } else {
        // Use corrected or original scaffolds (depending on Inspector setting)
        ch_final_scaffolds = ch_scaffolds_for_decontam
    }
    
    /*
    ========================================================================================
        STEP 10-12: Conditional Second Round of Scaffolding
        Only runs if:
        - Inspector correction on scaffolds is enabled, OR
        - Decontamination on scaffolds is enabled, OR
        - Explicitly enabled via --run_scaffold_round2 true
        
        Can be explicitly disabled via --run_scaffold_round2 false
    ========================================================================================
    */
    // Defaulted BEFORE the branch: assigned only inside it, this would not be
    // visible to CHIMERA_JOINS further down.
    ch_scaffold_round2_agp = Channel.empty()

    if (params.run_scaffold_round2) {
        
        log.info "[INFO] Running second round of scaffolding (scaffold correction or decontamination was performed)"
        
        /*
        ================================================================================
            STEP 10: Map Hi-C to Final Scaffolds
        ================================================================================
        */
        ch_final_scaffolds
            .map { meta, fasta -> [ meta.sample, meta, fasta ] }
            .combine( ch_hic_reads.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
            .map { sample, meta, fasta, hic_r1, hic_r2 ->
                tuple(meta, fasta, hic_r1, hic_r2, "scaffold")
            }
            .set { ch_hic_scaffold_mapping_input }

        MAP_HIC_TO_SCAFFOLD(ch_hic_scaffold_mapping_input)
        ch_versions = ch_versions.mix(MAP_HIC_TO_SCAFFOLD.out.versions)

        // checkpoint: scaffold_round2_raw_map (BAM-level only)
        MAP_HIC_TO_SCAFFOLD.out.bam
            .map { meta, stage, bam, bai -> tuple(meta, "scaffold_round2_raw_map", bam, bai) }
            .set { ch_hic_scaffold_raw_bam_for_qc }

        HIC_BAM_METRICS_SCAFFOLD(ch_hic_scaffold_raw_bam_for_qc)
        ch_all_bam_metrics = ch_all_bam_metrics.mix(HIC_BAM_METRICS_SCAFFOLD.out.metrics)

        /*
        ================================================================================
            STEP 11: Filter Hi-C BAM mapped to final scaffolds
        ================================================================================
        */
        MAP_HIC_TO_SCAFFOLD.out.bam
            .join(ch_final_scaffolds)
            .map { meta, stage, bam, bai, scaffold_fasta ->
                tuple(meta, stage, bam, bai, scaffold_fasta)
            }
            .set { ch_bam_with_scaffold }

        FILTER_HIC_BAM_SCAFFOLD(ch_bam_with_scaffold)
        ch_versions = ch_versions.mix(FILTER_HIC_BAM_SCAFFOLD.out.versions)

        // checkpoint: scaffold_round2_filtered (pairs-level + retention); already in scaffold1 names -> []
        FILTER_HIC_BAM_SCAFFOLD.out.pairs
            .join(FILTER_HIC_BAM_SCAFFOLD.out.parse_stats)
            .join(FILTER_HIC_BAM_SCAFFOLD.out.dedup_stats)
            .map { meta, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats ->
                tuple(meta, "scaffold_round2_filtered", pairs_gz, [], parse_stats, dedup_stats)
            }
            .set { ch_hic_pairs_scaffold_round2_filtered_for_qc }

        HIC_PAIRS_METRICS_SCAFFOLD(ch_hic_pairs_scaffold_round2_filtered_for_qc)
        ch_all_pairs_metrics = ch_all_pairs_metrics.mix(HIC_PAIRS_METRICS_SCAFFOLD.out.metrics)


        /*
        ================================================================================
            STEP 12: Second Round of Scaffolding
        ================================================================================
        */
        FILTER_HIC_BAM_SCAFFOLD.out.bam
            .join(ch_final_scaffolds)
            .map { meta, stage, bam, bai, scaffold_fasta ->
                tuple(meta, bam, bai, scaffold_fasta, "round2", [ min_contig_length: params.yahs_round2_min_contig_bp, min_mapq: params.yahs_round2_min_mapq, resolutions: params.yahs_round2_resolutions, rounds_per_resolution: params.yahs_round2_rounds_per_resolution, enzyme: params.yahs_round2_enzyme, no_contig_ec: !params.yahs_round2_contig_ec, no_scaffold_ec: !params.yahs_round2_scaffold_ec ])
            }
            .set { ch_second_scaffolding_input }

        SCAFFOLD_HIC_ROUND2(ch_second_scaffolding_input)
        ch_final_scaffolds_round2 = SCAFFOLD_HIC_ROUND2.out.scaffolds   // per-hap (meta, fasta)
        ch_scaffold_round2_agp    = SCAFFOLD_HIC_ROUND2.out.agp

        // checkpoint: scaffold_round2_space (relabel scaffold1->scaffold2 via round2 AGP; no remap)
        FILTER_HIC_BAM_SCAFFOLD.out.pairs
            .join(SCAFFOLD_HIC_ROUND2.out.agp)
            .join(FILTER_HIC_BAM_SCAFFOLD.out.parse_stats)
            .join(FILTER_HIC_BAM_SCAFFOLD.out.dedup_stats)
            .map { meta, stage, pairs_gz, agp, stage2, parse_stats, stage3, dedup_stats ->
                tuple(meta, "scaffold_round2_space", pairs_gz, agp, parse_stats, dedup_stats)
            }
            .set { ch_hic_pairs_scaffold_round2_space_for_qc }

        HIC_PAIRS_METRICS_SCAFFOLDSCAF(ch_hic_pairs_scaffold_round2_space_for_qc)
        ch_all_pairs_metrics = ch_all_pairs_metrics.mix(HIC_PAIRS_METRICS_SCAFFOLDSCAF.out.metrics)
        
    } else {
        log.info "[INFO] Skipping second round of scaffolding (no scaffold correction or decontamination)"
        
        // Create empty channel for round 2 scaffolds when not running
        ch_final_scaffolds_round2 = Channel.empty()
    }

    /*
    ========================================================================================
        STEP 13: Gap Filling
        Fills gaps in final scaffolded assemblies using HiFi reads
        Operates on the final scaffold output from either:
        - Round 2 scaffolding (if round 2 was run)
        - Decontaminated scaffolds (if decontamination on scaffolds was run)
        - Corrected scaffolds (if correction on scaffolds was run)
        - Original scaffolds (from round 1)
    ========================================================================================
    */
    
    // Determine which scaffolds to gap fill based on what was run
    if (params.run_scaffold_round2) {
        // Use round 2 scaffolds
        ch_scaffolds_for_gap_filling = ch_final_scaffolds_round2
    } else {
        // Use round 1 final scaffolds (corrected/decontaminated if those options were chosen)
        ch_scaffolds_for_gap_filling = ch_final_scaffolds
    }

    // HiFi-only assemblies are NOT gap-filled — no Hi-C scaffolding means no scaffold gaps to
    // close. They rejoin the finishing chain at teloclip/finalize below (they still have HiFi
    // reads), mirroring how short-read rejoins at ch_final_assembly.
    
    // Combine scaffolds with sample HiFi reads for gap filling (key on meta.sample)
    ch_scaffolds_for_gap_filling
        .map { meta, scaffold -> [ meta.sample, meta, scaffold ] }
        .combine( ch_hifi_reads.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
        .map { sample, meta, scaffold, hifi_fastq -> tuple(meta, scaffold, hifi_fastq) }
        .set { ch_gap_filling_input }
    
    // Run gap filling
    GAP_FILLING(ch_gap_filling_input)
    ch_versions = ch_versions.mix(GAP_FILLING.out.versions)

    // Gap-filled Hi-C scaffolds + HiFi-only contigs (which correctly skipped gap-fill) both
    // continue to teloclip/finalize.
    ch_round1 = SCAFFOLD_HIC_ROUND1.out.scaffolds
    ch_round1_agp = SCAFFOLD_HIC_ROUND1.out.agp
    ch_contig_pairs = FILTER_HIC_BAM.out.pairs
    ch_filled = GAP_FILLING.out.filled_assembly
    }
    ch_post_gap_fill = ch_filled.mix(ch_hifi_only_scaffolds)

    /*
    ========================================================================================
        STEP 13b: Teloclip — Extend scaffolds with missing telomeres (Optional)
        Maps raw HiFi reads back to gap-filled scaffolds to find soft-clipped
        alignments at scaffold ends containing telomeric motifs, then appends
        the overhang sequence to recover missing telomeres.
    ========================================================================================
    */
    if (capabilities.hifiasm && params.run_teloclip_extend) {
        // Combine gap-filled assemblies with sample HiFi reads (key on meta.sample)
        ch_post_gap_fill
            .map { meta, filled_fa -> [ meta.sample, meta, filled_fa ] }
            .combine( ch_hifi_reads.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
            .map { sample, meta, filled_fa, hifi_fastq -> tuple(meta.taxid?.toString(), meta, filled_fa, hifi_fastq) }
            .combine( ch_telo_by_taxid, by: 0 )
            .map { taxid, meta, filled_fa, hifi_fastq, telo -> tuple(meta, filled_fa, hifi_fastq, telo) }
            .set { ch_teloclip_input }

        TELOCLIP_EXTEND(ch_teloclip_input)
        ch_versions = ch_versions.mix(TELOCLIP_EXTEND.out.versions)

        // Collect teloclip stats across all haplotypes
        COLLECT_TELOCLIP_STATS(
            TELOCLIP_EXTEND.out.stats.map { meta, stats -> stats }.collect()
        )

        // The teloclip-extended assembly becomes the "final" assembly
        ch_extended = TELOCLIP_EXTEND.out.extended_assembly
        ch_final_assembly = ch_extended
        ch_teloclip_stats_for_report = COLLECT_TELOCLIP_STATS.out.stats.ifEmpty(file('NO_TELOCLIP'))
    } else {
        // No teloclip — gap-filled assembly IS the final assembly
        ch_final_assembly = ch_post_gap_fill
        ch_teloclip_stats_for_report = Channel.of(file('NO_TELOCLIP'))
    }


    emit:
    assembly = ch_final_assembly
    shortread = ch_shortread_finished
    round1 = ch_round1
    round1_agp = ch_round1_agp
    round2 = ch_final_scaffolds_round2
    round2_agp = ch_scaffold_round2_agp
    contig_pairs = ch_contig_pairs
    corrected = ch_scaffold_corrected
    decontaminated = ch_scaffold_decontaminated
    filled = ch_filled
    extended = ch_extended
    bam_metrics = ch_all_bam_metrics
    pairs_metrics = ch_all_pairs_metrics
    teloclip_stats = ch_teloclip_stats_for_report
    versions = ch_versions
}
