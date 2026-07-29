/*
========================================================================================
    FINAL_HIC_MAPS  —  Hi-C mapping + contact maps (+ optional compartments/TADs)
                       on the finalized assembly
========================================================================================
    Repo location:  workflows/final_hic_maps.nf

    Extracted verbatim from main.nf STEP 14 (the `if (params.run_final_contact_maps)`
    block). The caller wraps the invocation in that same guard, so this subworkflow
    assumes it should run; the internal `params.run_hic_balance` guard still controls
    whether compartments/TADs are computed.

    Body is unchanged except the trimmed-Hi-C input, which was `TRIM_HIC.out.trimmed_reads`
    in main and is now the `ch_trimmed_hic` take: channel.

    NOTE: matches current behaviour by NOT emitting a versions channel — the final
    mapping/contact processes were not mixed into ch_versions in main. Add a `versions`
    emit here when wiring software-version capture (piece 5).
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
    ch_trimmed_hic            //  TRIM_HIC.out.trimmed_reads        tuple(meta, r1, r2)
    ch_compartments_script    //  file: py_scripts/plot_compartments_pc1_genomewide.py
    ch_tad_book_script        //  file: py_scripts/make_tad_book.py

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

    FILTER_HIC_BAM_FINAL.out.pairs
        .join(ch_finalized_assembly)
        .map { meta, stage, pairs_gz, final_fasta ->
            tuple(meta, pairs_gz, final_fasta, "final")
        }
        .set { ch_contact_map_final_input }

    CONTACT_MAP_FINAL(ch_contact_map_final_input)

    if (params.run_hic_balance) {
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

    emit:
    contact_maps  = CONTACT_MAP_FINAL.out.contact_maps
    bam_metrics   = HIC_BAM_METRICS_FINAL.out.metrics
    pairs_metrics = HIC_PAIRS_METRICS_FINAL.out.metrics
}
