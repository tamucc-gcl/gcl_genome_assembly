include { CACTUS_PANGENOME } from '../modules/cactus_pangenome.nf'
include { PANGENOME_PREFLIGHT } from '../modules/pangenome_preflight.nf'
include { PANGENOME_GRAPH_CONTRACT } from '../modules/pangenome_graph_contract.nf'
include { PANGENOME_INPUT_AUDIT } from '../modules/pangenome_input_audit.nf'

workflow PANGENOME_CONSTRUCTION {
    take:
    ch_cohorts  // tuple(taxid, reference_name, names, fastas, ref_contigs, cohort)
    main:
    ch_clip_stats_input = Channel.empty()
    ch_report = Channel.empty()
    ch_manifest = Channel.empty()
    ch_versions = Channel.empty()
    PANGENOME_PREFLIGHT(ch_cohorts.map { row -> true }.unique())
    PANGENOME_INPUT_AUDIT(ch_cohorts,
        file("${projectDir}/py_scripts/pangenome_input_audit.py", checkIfExists: true))
    ch_report = PANGENOME_INPUT_AUDIT.out.report
    if (!params.pangenome_validate_only) {
    ch_checked = ch_cohorts.join(PANGENOME_INPUT_AUDIT.out.ready)
        .map { taxid, ref, names, fastas, contigs, cohort, audit -> tuple(taxid, ref, names, fastas, contigs) }
    CACTUS_PANGENOME(ch_checked.combine(PANGENOME_PREFLIGHT.out.ready.first())
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
    ch_report = PANGENOME_GRAPH_CONTRACT.out.report
    ch_manifest = PANGENOME_GRAPH_CONTRACT.out.manifest
    ch_versions = CACTUS_PANGENOME.out.versions
    }
    emit:
    clip_stats_input = ch_clip_stats_input
    report = ch_report
    manifest = ch_manifest
    identities = PANGENOME_INPUT_AUDIT.out.identities
    versions = ch_versions
}
