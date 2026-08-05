/*
 * FINAL_VIZ — cross-assembly comparison + telomere visualization for the finalized
 * per-hap assemblies. Extracted from main.nf to keep the entry workflow under Groovy's
 * 64 KB compile limit. Pure viz/QC: assemblies in, plots + summaries out.
 */
include { SETUP_PAFR; PAIRWISE_ALIGNMENT; COLLECT_PAIRWISE_RESULTS } from '../modules/pairwise_alignment.nf'
include { RIPARIAN_PLOT }                                            from '../modules/riparian_plot.nf'
include { QUAST_FINAL }                                     from '../modules/quast.nf'
include { TIDK_EXPLORE; TIDK_SEARCH; TIDK_PLOT; TIDK_SUMMARIZE; COLLECT_TIDK_RESULTS } from '../modules/tidk.nf'

workflow FINAL_VIZ {
    take:
    ch_final_by_id        // tuple(haplotype_id, fasta)
    ch_final_by_id_telo   // tuple(haplotype_id, fasta, telomere_motif)
    ch_dotplot_script     // path
    ch_riparian_script    // path

    main:
    // ---- Pairwise dotplots + riparian ----
    if (params.run_pairwise_alignments) {
        SETUP_PAFR()
        ch_final_by_id
            .toSortedList { a, b -> a[0] <=> b[0] }
            .flatMap { assemblies ->
                def pairs = []
                if (params.pairwise_alignment_mode == 'within_sample') {
                    def grouped = assemblies.groupBy { haplotype_id, fasta ->
                        haplotype_id.replaceAll(/_(hap[12]|primary)$/, '')
                    }
                    grouped.each { sample_id, haps ->
                        if (haps.size() == 2) {
                            def sorted = haps.sort { it[0] }
                            pairs << tuple(sorted[0][0], sorted[0][1], sorted[1][0], sorted[1][1])
                        }
                    }
                } else {
                    for (int i = 0; i < assemblies.size(); i++) {
                        for (int j = i + 1; j < assemblies.size(); j++) {
                            def (id1, fa1) = assemblies[i]
                            def (id2, fa2) = assemblies[j]
                            pairs << tuple(id1, fa1, id2, fa2)
                        }
                    }
                }
                return pairs
            }
            .set { ch_pairwise_input }

        PAIRWISE_ALIGNMENT(ch_pairwise_input, SETUP_PAFR.out.ready, ch_dotplot_script)

        ch_riparian_input = PAIRWISE_ALIGNMENT.out.paf
            .map { hap1, hap2, paf -> tuple(hap1, hap2, paf) }
            .combine(ch_final_by_id, by: 0)
            .map { hap1, hap2, paf, fasta1 -> tuple(hap2, hap1, fasta1, paf) }
            .combine(ch_final_by_id, by: 0)
            .map { hap2, hap1, fasta1, paf, fasta2 -> tuple(hap1, fasta1, hap2, fasta2, paf) }

        RIPARIAN_PLOT(ch_riparian_input, ch_riparian_script)

        COLLECT_PAIRWISE_RESULTS(
            PAIRWISE_ALIGNMENT.out.qc.map { id1, id2, qc_file -> qc_file }.collect()
        )
        ch_pw_summary  = COLLECT_PAIRWISE_RESULTS.out.summary.ifEmpty(file('NO_PAIRWISE'))
        ch_pw_dotplot  = PAIRWISE_ALIGNMENT.out.dotplot
        ch_pw_riparian = RIPARIAN_PLOT.out.riparian
    } else {
        ch_pw_summary  = Channel.of(file('NO_PAIRWISE'))
        ch_pw_dotplot  = Channel.empty()
        ch_pw_riparian = Channel.empty()
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
    QUAST_FINAL(ch_quast_assemblies.collect(), ch_quast_labels.flatten().collect())

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
    versions         = QUAST_FINAL.out.versions.mix(TIDK_EXPLORE.out.versions)
}