/*
========================================================================================
    FINAL_HIC_MAPS  —  Hi-C mapping + contact maps (+ optional compartments/TADs)
                       on the finalized assembly
========================================================================================
    Repo location:  workflows/final_hic_maps.nf

    Core final QC: map and filter once, then share filtered pairs between metrics
    and contact maps. The caller enables this for final_only/all_stages QC.
    run_final_contact_maps controls plots, not evidence generation. Compartments
    and TADs additionally require post-assembly analyses to be enabled.
========================================================================================
*/

include { MAP_HIC_TO_ASSEMBLY as MAP_HIC_TO_FINAL } from '../modules/map_hic_to_assembly.nf'
include { FILTER_HIC_BAM as FILTER_HIC_BAM_FINAL } from '../modules/filter_hic_bam.nf'
include { HIC_BAM_METRICS as HIC_BAM_METRICS_FINAL; HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_FINAL } from '../modules/hic_mapping_metrics.nf'
include { CONTACT_MAP as CONTACT_MAP_FINAL } from '../modules/contact_map.nf'
include { HIC_COMPARTMENTS } from '../modules/hic_compartments.nf'
include { HIC_TADS } from '../modules/hic_tads.nf'

workflow FINAL_HIC_MAPS {

    take:
    ch_finalized_assembly     //  FINALIZE_ASSEMBLY.out.assembly   tuple(meta, fasta)
    ch_trimmed_hic            // per sample, ordered mate lists and library metadata
    ch_compartments_script    //  file: py_scripts/plot_compartments_pc1_genomewide.py
    ch_tad_book_script        //  file: py_scripts/make_tad_book.py
    post_analyses_on          // optional compartments/TADs remain post-assembly analyses

    main:
    // Combine final per-hap assemblies with sample Hi-C reads (key on meta.sample)
    ch_finalized_assembly
        .map { meta, final_fa -> [ meta.sample, meta, final_fa ] }
        .combine( ch_trimmed_hic.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
        .map { sample, meta, final_fa, hic_r1, hic_r2 ->
            tuple(meta, final_fa, hic_r1, hic_r2, 'final')
        }
        .set { ch_final_hic_map_input }

    MAP_HIC_TO_FINAL(ch_final_hic_map_input)

    // checkpoint: final_raw_map (BAM-level only)
    MAP_HIC_TO_FINAL.out.bam
        .map { meta, stage, bam, bai -> tuple(meta, "final_raw_map", bam, bai) }
        .set { ch_final_raw_bam_qc }

    HIC_BAM_METRICS_FINAL(ch_final_raw_bam_qc)

    MAP_HIC_TO_FINAL.out.bam
        .join(ch_finalized_assembly)
        .map { meta, stage, bam, bai, final_fa ->
            tuple(meta, "final", bam, bai, final_fa)
        }
        .set { ch_final_filter_input }

    FILTER_HIC_BAM_FINAL(ch_final_filter_input)

    // checkpoint: final_filtered (pairs-level + retention); already in final names -> []
    FILTER_HIC_BAM_FINAL.out.pairs
        .join(FILTER_HIC_BAM_FINAL.out.parse_stats)
        .join(FILTER_HIC_BAM_FINAL.out.dedup_stats)
        .map { meta, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats ->
            tuple(meta, "final_filtered", pairs_gz, [], parse_stats, dedup_stats)
        }
        .set { ch_final_pairs_qc }

    HIC_PAIRS_METRICS_FINAL(ch_final_pairs_qc)

    ch_maps = Channel.empty()
    ch_versions = MAP_HIC_TO_FINAL.out.versions.mix(FILTER_HIC_BAM_FINAL.out.versions)
    if (params.run_final_contact_maps != false) {
        FILTER_HIC_BAM_FINAL.out.pairs
            .join(ch_finalized_assembly)
            .map { meta, stage, pairs_gz, final_fasta ->
                tuple(meta, pairs_gz, final_fasta, "final")
            }
            .set { ch_contact_map_final_input }

        CONTACT_MAP_FINAL(ch_contact_map_final_input)
        ch_maps = CONTACT_MAP_FINAL.out.contact_maps
        ch_versions = ch_versions.mix(CONTACT_MAP_FINAL.out.versions)

        if (post_analyses_on && params.run_hic_balance) {
            HIC_COMPARTMENTS(
                CONTACT_MAP_FINAL.out.mcool,
                params.compartment_resolution ?: 250000,
                params.compartment_min_contig_bp ?: 5000000,
                params.compartment_max_contigs ?: 30,
                ch_compartments_script
            )

            HIC_TADS(
                CONTACT_MAP_FINAL.out.mcool,
                params.tad_resolution ?: 50000,
                params.tad_window_bp ?: 500000,
                params.tad_min_contig_bp ?: 5000000,
                params.tad_max_contigs ?: 0,
                ch_tad_book_script
            )
        }
    }

    emit:
    contact_maps  = ch_maps
    bam_metrics   = HIC_BAM_METRICS_FINAL.out.metrics
    pairs_metrics = HIC_PAIRS_METRICS_FINAL.out.metrics
    bam = MAP_HIC_TO_FINAL.out.bam
    filtered_bam = FILTER_HIC_BAM_FINAL.out.bam
    pairs = FILTER_HIC_BAM_FINAL.out.pairs
    versions = ch_versions
}
