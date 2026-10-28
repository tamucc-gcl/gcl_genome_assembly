/*
 * FINAL_VIZ — cross-assembly comparison + telomere visualization for the finalized
 * per-hap assemblies. Extracted from main.nf to keep the entry workflow under Groovy's
 * 64 KB compile limit. Pure viz/QC: assemblies in, plots + summaries out.
 */
include { PAIRWISE_COMPARISONS } from './pairwise_comparisons.nf'
include { QUAST_FINAL }                                     from '../modules/quast.nf'
include { TIDK_EXPLORE; TIDK_SEARCH; TIDK_PLOT; TIDK_SUMMARIZE; COLLECT_TIDK_RESULTS } from '../modules/tidk.nf'

workflow FINAL_VIZ {
    take:
    ch_finalized         // tuple(meta, final FASTA, final FAI)
    ch_final_by_id_telo   // tuple(haplotype_id, fasta, telomere_motif)

    main:
    ch_final_by_id=ch_finalized.map { meta,fa,fai -> tuple(meta.id,fa) }
    ch_pairwise_versions=Channel.empty()
    if (params.run_pairwise_alignments) {
        PAIRWISE_COMPARISONS(ch_finalized)
        ch_pw_summary=PAIRWISE_COMPARISONS.out.summary
        ch_pw_dotplot=PAIRWISE_COMPARISONS.out.dotplot
        ch_pw_riparian=PAIRWISE_COMPARISONS.out.riparian
        ch_pairwise_versions=PAIRWISE_COMPARISONS.out.versions
    } else {
        ch_pw_summary=Channel.value(file("${projectDir}/assets/NO_PAIRWISE",checkIfExists:true))
        ch_pw_dotplot=Channel.empty()
        ch_pw_riparian=Channel.empty()
    }

    // ---- Cross-sample QUAST comparison ----
    ch_final_by_id
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { list ->
            def labels     = list.collect { it[0] }
            def assemblies = list.collect { it[1] }
            tuple(assemblies, labels)
        }
        .set { ch_quast_final_input }
    ch_quast_final_input.map { assemblies, labels -> assemblies }.set { ch_quast_assemblies }
    ch_quast_final_input.map { assemblies, labels -> labels }.set { ch_quast_labels }
    QUAST_FINAL(ch_quast_assemblies, ch_quast_labels)

    // ---- Telomere detection ----
    TIDK_EXPLORE(ch_final_by_id)
    TIDK_SEARCH(ch_final_by_id_telo)
    TIDK_PLOT(TIDK_SEARCH.out.search_tsv)
    TIDK_SUMMARIZE(TIDK_SEARCH.out.search_tsv)
    COLLECT_TIDK_RESULTS(
        TIDK_SUMMARIZE.out.summary.map   { haplotype_id, summary   -> summary   }.collect(),
        TIDK_SUMMARIZE.out.telomeres.map { haplotype_id, telomeres -> telomeres }.collect()
    )

    emit:
    pairwise_summary = ch_pw_summary
    dotplot          = ch_pw_dotplot
    riparian         = ch_pw_riparian
    tidk_plot        = TIDK_PLOT.out.plot
    telomere_summary = COLLECT_TIDK_RESULTS.out.summary
    versions         = QUAST_FINAL.out.versions.mix(TIDK_EXPLORE.out.versions).mix(ch_pairwise_versions)
}