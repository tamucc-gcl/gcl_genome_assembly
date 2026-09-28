include { PANGENOME_STATS } from '../modules/pangenome_stats.nf'

/* Published-tool sharing/ordination enters in the next batch.
 * Only validated CLIP handles enter the biological statistics workflow.
 */
workflow PANGENOME_ANALYSIS {
    take:
    ch_clip_stats_input
    main:
    PANGENOME_STATS(ch_clip_stats_input)
    emit:
    stats = PANGENOME_STATS.out.vg_stats
    versions = Channel.empty()
}