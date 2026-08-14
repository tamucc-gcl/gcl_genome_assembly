/*
========================================================================================
    HIC_PREFINAL_PAIRS  —  map + filter Hi-C on the PRE-finalize assembly
========================================================================================
    Repo location: workflows/hic_prefinal_pairs.nf

    Extracted from FINAL_HIC_MAPS and moved UPSTREAM of HARMONIZE_SCAFFOLDS, so the pairs
    exist in the same coordinate system as the assemblies being harmonized. That is what
    makes the Hi-C junction evidence channel possible; nothing else changed about the
    mapping itself.

    COORDINATES AND PUBLISHING
      BAMs stay in pre-finalize coordinates and are published as such, under
      bam/hic/prefinal/. A BAM cannot be lifted the way pairs can -- it carries
      reverse-complemented sequence and flags, not just coordinates -- so there is no
      finalized BAM anywhere in the pipeline. The pairs ARE lifted, by
      HIC_LIFT_HARMONIZED_PAIRS, and published in finalized coordinates under
      bam/hic/final/filtered/. Prefinal pairs are published alongside the prefinal BAMs;
      the directory name is the coordinate system.

    One Hi-C alignment total, unchanged.

    Take : ch_assembly    tuple(meta, fasta)   post-teloclip (ch_final_assembly)
           ch_trimmed_hic tuple(meta, r1, r2)
    Emit : pairs / bam / bam_metrics / pairs_metrics / versions
========================================================================================
*/

include { MAP_HIC_TO_ASSEMBLY as MAP_HIC_PREFINAL } from '../modules/map_hic_to_assembly.nf'
include { FILTER_HIC_BAM as FILTER_HIC_BAM_PREFINAL } from '../modules/filter_hic_bam.nf'
include { HIC_BAM_METRICS as HIC_BAM_METRICS_PREFINAL; HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_PREFINAL } from '../modules/hic_mapping_metrics.nf'

workflow HIC_PREFINAL_PAIRS {

    take:
    ch_assembly
    ch_trimmed_hic

    main:
    ch_versions = Channel.empty()

    ch_assembly
        .map { meta, fa -> [ meta.sample, meta, fa ] }
        .combine( ch_trimmed_hic.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
        .map { sample, meta, fa, r1, r2 -> tuple(meta, fa, r1, r2, 'prefinal') }
        .set { ch_map_input }

    MAP_HIC_PREFINAL(ch_map_input)
    ch_versions = ch_versions.mix(MAP_HIC_PREFINAL.out.versions)

    MAP_HIC_PREFINAL.out.bam
        .map { meta, stage, bam, bai -> tuple(meta, "final_raw_map", bam, bai) }
        .set { ch_bam_qc }
    HIC_BAM_METRICS_PREFINAL(ch_bam_qc)

    MAP_HIC_PREFINAL.out.bam
        .join(ch_assembly)
        .map { meta, stage, bam, bai, fa -> tuple(meta, 'prefinal', bam, bai, fa) }
        .set { ch_filter_input }

    FILTER_HIC_BAM_PREFINAL(ch_filter_input)
    ch_versions = ch_versions.mix(FILTER_HIC_BAM_PREFINAL.out.versions)

    FILTER_HIC_BAM_PREFINAL.out.pairs
        .join(FILTER_HIC_BAM_PREFINAL.out.parse_stats)
        .join(FILTER_HIC_BAM_PREFINAL.out.dedup_stats)
        .map { meta, s1, pairs_gz, s2, parse_stats, s3, dedup_stats ->
            tuple(meta, "final_filtered", pairs_gz, [], parse_stats, dedup_stats)
        }
        .set { ch_pairs_qc }
    HIC_PAIRS_METRICS_PREFINAL(ch_pairs_qc)

    emit:
    pairs         = FILTER_HIC_BAM_PREFINAL.out.pairs      // tuple(meta, stage, pairs_gz)
    bam           = MAP_HIC_PREFINAL.out.bam
    bam_metrics   = HIC_BAM_METRICS_PREFINAL.out.metrics
    pairs_metrics = HIC_PAIRS_METRICS_PREFINAL.out.metrics
    versions      = ch_versions
}
