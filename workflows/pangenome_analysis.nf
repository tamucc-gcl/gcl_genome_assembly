include { PANGENOME_STATS } from '../modules/pangenome_stats.nf'
include { PANGENOME_SHARING } from '../modules/pangenome_sharing.nf'
include { PANGENOME_SHARING_PLOTS } from '../modules/pangenome_sharing_plots.nf'

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
    }
    emit:
    stats = PANGENOME_STATS.out.vg_stats
    versions = ch_versions
    report = ch_report
}
