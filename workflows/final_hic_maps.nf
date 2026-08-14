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
include { HIC_LIFT_HARMONIZED_PAIRS } from '../modules/hic_lift_harmonized_pairs.nf'
include { HIC_COMPARTMENTS } from '../modules/hic_compartments.nf'
include { HIC_TADS } from '../modules/hic_tads.nf'

workflow FINAL_HIC_MAPS {

    take:
    ch_finalized_assembly     //  FINALIZE_ASSEMBLY.out.assembly       tuple(meta, fasta)
    ch_prefinalize_assembly   //  ch_final_assembly, post-teloclip     tuple(meta, fasta)
    ch_applied_lift           //  FINALIZE_ASSEMBLY.out.applied_lift   tuple(meta, tsv)
    ch_trimmed_hic            //  TRIM_HIC.out.trimmed_reads           tuple(meta, r1, r2)
    ch_compartments_script    //  file: py_scripts/plot_compartments_pc1_genomewide.py
    ch_tad_book_script        //  file: py_scripts/make_tad_book.py

    main:
    // Declared here rather than in main.nf: main.nf sits ~1 KB below Groovy's 65,535-char
    // compiled-unit limit, so lines there are the scarce resource.
    ch_hic_lift_script = file("${projectDir}/py_scripts/lift_harmonized_pairs.py",
                              checkIfExists: true)

    // Map Hi-C to the PRE-finalize (post-teloclip) assembly, then lift the pairs into
    // finalized coordinates arithmetically. Same single mapping as before -- but producing
    // the pairs before renaming is what lets HARMONIZE consume Hi-C contacts without a
    // dependency cycle, which is the whole point of the reorder.
    ch_prefinalize_assembly
        .map { meta, pre_fa -> [ meta.sample, meta, pre_fa ] }
        .combine( ch_trimmed_hic.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
        .map { sample, meta, pre_fa, hic_r1, hic_r2 ->
            tuple(meta, pre_fa, hic_r1, hic_r2, 'final')
        }
        .set { ch_final_hic_map_input }

    MAP_HIC_TO_FINAL(ch_final_hic_map_input)

    // checkpoint: final_raw_map (BAM-level only)
    MAP_HIC_TO_FINAL.out.bam
        .map { meta, stage, bam, bai -> tuple(meta, "final_raw_map", bam, bai) }
        .set { ch_final_raw_bam_qc }

    HIC_BAM_METRICS_FINAL(ch_final_raw_bam_qc)

    // the BAM is in PRE-finalize coordinates, so filtering must see that fasta
    MAP_HIC_TO_FINAL.out.bam
        .join(ch_prefinalize_assembly)
        .map { meta, stage, bam, bai, pre_fa ->
            tuple(meta, "final", bam, bai, pre_fa)
        }
        .set { ch_final_filter_input }

    FILTER_HIC_BAM_FINAL(ch_final_filter_input)

    // ---- lift pairs into finalized coordinates (rename + revcomp; no remapping) --------
    FILTER_HIC_BAM_FINAL.out.pairs
        .join(ch_applied_lift)
        .map { meta, stage, pairs_gz, lift_tsv ->
            tuple(meta, "final", pairs_gz, lift_tsv)
        }
        .set { ch_lift_input }

    HIC_LIFT_HARMONIZED_PAIRS(ch_lift_input, ch_hic_lift_script)
    ch_final_pairs = HIC_LIFT_HARMONIZED_PAIRS.out.pairs

    // checkpoint: final_filtered. Metrics run on the LIFTED pairs so any per-scaffold
    // breakdown reports finalized names; the aggregate retention numbers are unchanged
    // either way, since the lift preserves every pair.
    ch_final_pairs
        .join(FILTER_HIC_BAM_FINAL.out.parse_stats)
        .join(FILTER_HIC_BAM_FINAL.out.dedup_stats)
        .map { meta, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats ->
            tuple(meta, "final_filtered", pairs_gz, [], parse_stats, dedup_stats)
        }
        .set { ch_final_pairs_qc }

    HIC_PAIRS_METRICS_FINAL(ch_final_pairs_qc)

    ch_final_pairs
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
    final_pairs   = ch_final_pairs                          // finalized coordinates
    lift_stats    = HIC_LIFT_HARMONIZED_PAIRS.out.stats
    lift_versions = HIC_LIFT_HARMONIZED_PAIRS.out.versions   // not yet mixed in main.nf
}
