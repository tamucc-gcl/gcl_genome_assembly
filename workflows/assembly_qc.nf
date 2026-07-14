/*
========================================================================================
    ASSEMBLY QC SUBWORKFLOW (CACHE-OPTIMIZED, META-THREADED)
========================================================================================
    Repo location: workflows/assembly_qc.nf

    Comprehensive QC for phased genome assemblies.
      - QUAST:      per-sample (both haplotypes)
      - MERQURY:    per-sample (both haplotypes) — uses pre-built meryl database
      - BUSCO:      per-haplotype
      - MAPPING_QC: per-haplotype (map HiFi reads back to assembly)
      - COMBINE:    aggregate + visualize per sample

    KEY CHANGE (meta refactor):
      - `take: assemblies` is now the per-haplotype `(meta, fasta)` stream (NOT a
        re-paired triple). Callers just pass their per-hap channel — no groupTuple
        prep block at each call site.
      - Haplotypes are re-paired *here*, once, via groupKey(meta.sample, meta.n_hap),
        which replaces every `groupTuple(by:0, size:2)` (no size:2 blocking / haploid
        hang) and every `replaceAll(/_hap[12]$/,'')` (sample comes from meta.sample).
      - Per-sample caching preserved: each sample groups + emits independently once its
        meta.n_hap haplotypes arrive.

    PHASE 1 SCOPE: diploid path (n_hap == 2). QUAST/MERQURY take exactly two fastas.
    Haploid (n_hap == 1) QUAST/MERQURY handling is Phase 2.
========================================================================================
*/

nextflow.enable.dsl = 2

include { QUAST }               from '../modules/quast.nf'
include { MERQURY }             from '../modules/merqury.nf'
include { BUSCO }               from '../modules/busco.nf'
include { MAPPING_QC }          from '../modules/mapping_qc.nf'
include { COMBINE_ASSEMBLY_QC } from '../modules/combine_assembly_qc.nf'
combine_qc_script    = file("${projectDir}/r_scripts/combine_individual_assembly_qc.R", checkIfExists: true)

workflow ASSEMBLY_QC {
    take:
    assemblies   // channel: tuple(meta, fasta)          — per-haplotype
    hifi_reads   // channel: tuple(meta, hifi_fastq)      — per-sample (sample-level meta)
    meryl_db     // channel: tuple(meta, meryl_db)        — per-sample (sample-level meta)
    busco_db     // value:   path to pre-downloaded BUSCO lineage database
    qc_label     // value:   label for output subfolder

    main:

    /*
    ========================================================================================
        Re-pair haplotypes per sample (groupKey — no size:2 block) for QUAST + MERQURY.
        Emits (sample_id, [fasta ordered by meta.haplotype]).
    ========================================================================================
    */
    assemblies
        .map { meta, fasta -> tuple(groupKey(meta.sample, meta.n_hap), meta, fasta) }
        .groupTuple()
        .map { key, metas, fastas ->
            def sample_id = metas[0].sample
            def ordered = [metas, fastas].transpose().sort { it[0].haplotype }
            tuple(sample_id,
                  ordered.collect { it[1] },                              // [fastas]  (1 haploid | 2 diploid)
                  ordered.collect { "${sample_id}.${it[0].haplotype}" })  // [labels]  sample.hap1/.hap2 | sample.primary
        }
        .set { ch_paired_assemblies }   // (sample_id, [fastas], [labels])

    // QUAST — per-sample; 1 (haploid) or 2 (diploid) assemblies + matching labels
    QUAST(ch_paired_assemblies)

    // MERQURY — join re-paired assemblies with the sample's meryl DB (by sample_id)
    ch_paired_assemblies
        .map { sample_id, fastas, labels -> tuple(sample_id, fastas) }
        .join( meryl_db.map { meta, db -> tuple(meta.sample, db) } )
        .set { ch_merqury_input }   // (sample_id, [fastas], db)

    MERQURY(ch_merqury_input)

    /*
    ========================================================================================
        BUSCO — per-haplotype (carries meta straight through)
    ========================================================================================
    */
    BUSCO(assemblies, busco_db)

    /*
    ========================================================================================
        MAPPING QC — per-haplotype; combine each hap with its sample's HiFi reads
        (key on meta.sample — hifi_reads is per-sample)
    ========================================================================================
    */
    assemblies
        .map { meta, fasta -> [ meta.sample, meta, fasta ] }
        .combine( hifi_reads.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
        .map { sample, meta, fasta, hifi_fastq -> tuple(meta, fasta, hifi_fastq) }
        .set { ch_mapping_input }

    MAPPING_QC(ch_mapping_input)

    /*
    ========================================================================================
        Regroup per-haplotype BUSCO + MAPPING results back to per sample (groupKey).
        Emits (sample_id, [haplotype_ids ordered], [results ordered]) — the ordered
        haplotype-id list is what COMBINE_ASSEMBLY_QC writes into its manifest.
    ========================================================================================
    */
    BUSCO.out.results
        .map { meta, results -> tuple(groupKey(meta.sample, meta.n_hap), meta, results) }
        .groupTuple()
        .map { key, metas, results ->
            def sample_id = metas[0].sample
            def ordered = [metas, results].transpose().sort { it[0].haplotype }
            tuple(sample_id, ordered.collect { it[0].id }, ordered.collect { it[1] })
        }
        .set { ch_busco_by_sample }   // (sample_id, [hap_ids], [busco_results])

    MAPPING_QC.out.results
        .map { meta, results -> tuple(groupKey(meta.sample, meta.n_hap), meta, results) }
        .groupTuple()
        .map { key, metas, results ->
            def sample_id = metas[0].sample
            def ordered = [metas, results].transpose().sort { it[0].haplotype }
            tuple(sample_id, ordered.collect { it[0].id }, ordered.collect { it[1] })
        }
        .set { ch_mapping_by_sample }  // (sample_id, [hap_ids], [mapping_results])

    /*
    ========================================================================================
        Join all per-sample QC (by sample_id only) and attach qc_label directly.
        Each sample proceeds independently — no qc_label in the join key.
    ========================================================================================
    */
    QUAST.out.results
        .join(MERQURY.out.results, by: 0)
        .join(ch_busco_by_sample, by: 0)
        .join(ch_mapping_by_sample, by: 0)
        .map { sample_id, quast_results, merqury_results,
               hap_ids_busco, busco_results,
               hap_ids_mapping, mapping_results ->
            tuple(sample_id, qc_label,
                  quast_results, merqury_results,
                  hap_ids_busco, busco_results,
                  hap_ids_mapping, mapping_results)
        }
        .set { ch_all_qc_labeled }

    COMBINE_ASSEMBLY_QC(ch_all_qc_labeled, combine_qc_script)

    emit:
    assembly_summary = COMBINE_ASSEMBLY_QC.out.summary   // (sample_id, qc_label, summary_tsv)
    busco_results    = BUSCO.out.results                 // per-haplotype (meta, results)
}
