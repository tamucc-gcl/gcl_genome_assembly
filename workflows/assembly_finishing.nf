include { GAP_FILLING } from '../modules/gap_filling.nf'
include { TELOCLIP_EXTEND; COLLECT_TELOCLIP_STATS } from '../modules/teloclip.nf'

/* Sequence finishing follows chimera assessment/cuts. Name assignments travel separately:
 * neither gap filling nor telomere extension may rename, split, or reorder the assignment.
 * FINALIZE_ASSEMBLY checks sequence IDs and refreshes lengths before applying the map.
 */
workflow ASSEMBLY_FINISHING {
    take:
    ch_assemblies       // tuple(meta, fasta, name_map), including short-read pass-through
    ch_hifi
    ch_telo_by_taxid
    capabilities

    main:
    ch_versions = Channel.empty()
    ch_filled = Channel.empty()
    ch_extended = Channel.empty()
    ch_teloclip_stats = Channel.value(file('NO_TELOCLIP'))
    ch_maps = ch_assemblies.map { meta, fa, nm -> tuple(meta.id, nm) }
    ch_input = ch_assemblies.map { meta, fa, nm -> tuple(meta, fa) }
    ch_after_gap = ch_input
    if (capabilities.scaffold) {
        GAP_FILLING(ch_input.filter { meta, fa -> meta.assembler == 'hifiasm' && meta.hic }
            .map { meta, fa -> tuple(meta.sample, meta, fa) }
            .combine(ch_hifi.map { meta, fq -> tuple(meta.sample, fq) }, by: 0)
            .map { sample, meta, fa, fq -> tuple(meta, fa, fq) })
        ch_filled = GAP_FILLING.out.filled_assembly
        ch_versions = ch_versions.mix(GAP_FILLING.out.versions)
        ch_after_gap = ch_filled.mix(ch_input.filter { meta, fa -> !(meta.assembler == 'hifiasm' && meta.hic) })
    }
    ch_finished = ch_after_gap
    if (capabilities.hifiasm && params.run_teloclip_extend) {
        TELOCLIP_EXTEND(ch_after_gap.filter { meta, fa -> meta.assembler == 'hifiasm' }
            .map { meta, fa -> tuple(meta.sample, meta, fa) }
            .combine(ch_hifi.map { meta, fq -> tuple(meta.sample, fq) }, by: 0)
            .map { sample, meta, fa, fq -> tuple(meta.taxid.toString(), meta, fa, fq) }
            .combine(ch_telo_by_taxid.map { taxid, motif -> tuple(taxid.toString(), motif) }, by: 0)
            .map { taxid, meta, fa, fq, motif -> tuple(meta, fa, fq, motif) })
        ch_extended = TELOCLIP_EXTEND.out.extended_assembly
        ch_versions = ch_versions.mix(TELOCLIP_EXTEND.out.versions)
        COLLECT_TELOCLIP_STATS(TELOCLIP_EXTEND.out.stats.map { meta, f -> f }.collect())
        ch_teloclip_stats = COLLECT_TELOCLIP_STATS.out.stats.ifEmpty(file('NO_TELOCLIP'))
        ch_finished = ch_extended.mix(ch_after_gap.filter { meta, fa -> meta.assembler != 'hifiasm' })
    }
    ch_output = ch_finished.map { meta, fa -> tuple(meta.id, meta, fa) }
        .join(ch_maps, failOnDuplicate: true)
        .map { id, meta, fa, nm -> tuple(meta, fa, nm) }

    emit:
    assemblies = ch_output
    filled = ch_filled
    extended = ch_extended
    teloclip_stats = ch_teloclip_stats
    versions = ch_versions
}