include { PANGENOME_STATS } from '../modules/pangenome_stats.nf'
include { PANGENOME_SHARING } from '../modules/pangenome_sharing.nf'
include { PANGENOME_SHARING_PLOTS } from '../modules/pangenome_sharing_plots.nf'
include { PANGENOME_SIMILARITY; PANGENOME_COMPARISON_PLOTS } from '../modules/pangenome_similarity.nf'
include { PANGENOME_NODE_COVERAGE; PANGENOME_REGIONAL } from '../modules/pangenome_regional.nf'
include { PANGENOME_REGIONAL_PLOTS } from '../modules/pangenome_regional_plots.nf'

/* Only validated CLIP handles enter biological statistics.
 * Identity grouping is audited before Panacus; ordination follows separately.
 */
workflow PANGENOME_ANALYSIS {
    take:
    ch_clip_stats_input
    ch_clip_sharing_input
    main:
    ch_versions = Channel.empty()
    ch_report = Channel.empty()
    PANGENOME_STATS(ch_clip_stats_input)
    if (params.pangenome_growth) {
        PANGENOME_SHARING(ch_clip_sharing_input,
            file("${projectDir}/py_scripts/pangenome_sharing.py", checkIfExists: true))
        ch_versions = PANGENOME_SHARING.out.versions
        ch_report = PANGENOME_SHARING.out.report
        PANGENOME_SHARING_PLOTS(PANGENOME_SHARING.out.files,
            file("${projectDir}/r_scripts/pangenome_sharing_plots.R", checkIfExists: true))
        ch_versions = ch_versions.mix(PANGENOME_SHARING_PLOTS.out.versions)
        ch_report = ch_report.mix(PANGENOME_SHARING_PLOTS.out.report)
        ch_regional_input = ch_clip_sharing_input.join(PANGENOME_SHARING.out.files)
            .map { taxid, gfa, identities, products ->
                def files = products instanceof Collection ? products : [products]
                def groups = files.find { it.name == 'haplotype.groups.tsv' }
                def hist = files.find { it.name == 'haplotype.hist.tsv' }
                if (!groups || !hist) error 'Missing audited sharing inputs for regional attribution'
                tuple(taxid, gfa, groups, hist)
            }
        PANGENOME_NODE_COVERAGE(ch_regional_input.map { taxid, gfa, groups, hist -> tuple(taxid, gfa, groups) })
        PANGENOME_REGIONAL(ch_regional_input.join(PANGENOME_NODE_COVERAGE.out.coverage),
            file("${projectDir}/py_scripts/pangenome_regional.py", checkIfExists: true),
            file("${projectDir}/tests/test_pangenome_regional.py", checkIfExists: true))
        ch_versions = ch_versions.mix(PANGENOME_NODE_COVERAGE.out.versions)
        ch_report = ch_report.mix(PANGENOME_REGIONAL.out.report)
        PANGENOME_REGIONAL_PLOTS(PANGENOME_REGIONAL.out.sharing,
            file("${projectDir}/r_scripts/pangenome_regional_plots.R", checkIfExists: true))
        ch_report = ch_report.mix(PANGENOME_REGIONAL_PLOTS.out.report)
        ch_versions = ch_versions.mix(PANGENOME_REGIONAL_PLOTS.out.versions)
        if (params.pangenome_popstruct) {
            ch_similarity_input = ch_clip_sharing_input.join(PANGENOME_SHARING.out.files)
                .map { taxid, gfa, identities, products ->
                    def files = products instanceof Collection ? products : [products]
                    def groups = files.find { it.name == 'haplotype.groups.tsv' }
                    if (!groups) error 'Missing audited haplotype groups for similarity'
                    tuple(taxid, gfa, identities, groups)
                }
            PANGENOME_SIMILARITY(ch_similarity_input)
            PANGENOME_COMPARISON_PLOTS(PANGENOME_SIMILARITY.out.matrix,
                file("${projectDir}/r_scripts/pangenome_comparison_plots.R", checkIfExists: true))
            ch_versions = ch_versions.mix(PANGENOME_SIMILARITY.out.versions)
                .mix(PANGENOME_COMPARISON_PLOTS.out.versions)
            ch_report = ch_report.mix(PANGENOME_COMPARISON_PLOTS.out.report)
        }
    }
    emit:
    stats = PANGENOME_STATS.out.vg_stats
    versions = ch_versions
    report = ch_report
}
