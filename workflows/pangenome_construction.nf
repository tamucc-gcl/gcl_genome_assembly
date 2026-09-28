include { CACTUS_PANGENOME } from '../modules/cactus_pangenome.nf'
include { PANGENOME_PREFLIGHT } from '../modules/pangenome_preflight.nf'
include { PANGENOME_GRAPH_CONTRACT } from '../modules/pangenome_graph_contract.nf'

workflow PANGENOME_CONSTRUCTION {
    take:
    ch_cohorts  // tuple(taxid, reference_name, names, fastas, ref_contigs)
    main:
    PANGENOME_PREFLIGHT(ch_cohorts.map { row -> true }.unique())
    CACTUS_PANGENOME(ch_cohorts.combine(PANGENOME_PREFLIGHT.out.ready.first())
        .map { taxid, ref, names, fastas, contigs, check -> tuple(taxid, ref, names, fastas, contigs, check) })
    ch_contract = CACTUS_PANGENOME.out.all.map { taxid, products ->
        def files = products instanceof Collection ? products.toList() : [products]
        def suffixes = [
            biological_gbz: '.gbz', biological_gfa: '.gfa.gz', biological_odgi: '.og',
            haplotype_index: '.hapl', snarls: '.snarls',
            variant_gbz: '.gref.gbz', variant_gfa: '.gref.gfa.gz',
            variants_standard: '.gref.vcf.gz', variants_standard_index: '.gref.vcf.gz.tbi',
            variants_raw: '.gref.raw.vcf.gz', variant_segments: '.gref.gref-segs.tsv.gz',
            full_qc_gfa: '.full.gfa.gz', clipping_stats: '.stats.tgz']
        def selected = suffixes.collectEntries { role, suffix ->
            [(role): files.find { it.name == taxid + suffix }] }
        tuple(taxid, selected.collectEntries { role, f -> [(role): f?.name ?: ''] },
              selected.values().findAll { it != null }.unique())
    }
    PANGENOME_GRAPH_CONTRACT(ch_contract)
    ch_clip_stats_input = CACTUS_PANGENOME.out.gbz.join(CACTUS_PANGENOME.out.og)
        .join(PANGENOME_GRAPH_CONTRACT.out.ready)
        .map { taxid, gbz, og, check -> tuple(taxid, gbz, og) }
    emit:
    clip_stats_input = ch_clip_stats_input
    report = PANGENOME_GRAPH_CONTRACT.out.report
    manifest = PANGENOME_GRAPH_CONTRACT.out.manifest
    versions = CACTUS_PANGENOME.out.versions
}