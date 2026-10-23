include { CACTUS_PANGENOME } from '../modules/cactus_pangenome.nf'
include { PANGENOME_PREFLIGHT } from '../modules/pangenome_preflight.nf'
include { PANGENOME_GRAPH_CONTRACT } from '../modules/pangenome_graph_contract.nf'
include { PANGENOME_INPUT_AUDIT } from '../modules/pangenome_input_audit.nf'

workflow PANGENOME_CONSTRUCTION {
    take:
    ch_cohorts  // tuple(taxid, reference_name, names, fastas, ref_contigs, cohort)
    main:
    ch_clip_stats_input = Channel.empty()
    ch_clip_sharing_input = Channel.empty()
    ch_variant_input = Channel.empty()
    ch_report = Channel.empty()
    ch_manifest = Channel.empty()
    ch_versions = Channel.empty()
    ch_support_qc = Channel.empty()
    ch_source_evidence = Channel.empty()
    ch_assignment_input = Channel.empty()
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
    PANGENOME_GRAPH_CONTRACT(ch_contract,
        file("${projectDir}/py_scripts/pangenome_output_audit.py", checkIfExists: true))
    ch_clip_stats_input = CACTUS_PANGENOME.out.gbz.join(CACTUS_PANGENOME.out.og)
        .join(PANGENOME_GRAPH_CONTRACT.out.ready)
        .map { taxid, gbz, og, check -> tuple(taxid, gbz, og) }
    ch_report = ch_report.mix(PANGENOME_GRAPH_CONTRACT.out.report)
    ch_clip_sharing_input = CACTUS_PANGENOME.out.gfa
        .join(PANGENOME_INPUT_AUDIT.out.identities)
        .join(PANGENOME_GRAPH_CONTRACT.out.ready)
        .map { taxid, gfa, identities, check -> tuple(taxid, gfa, identities) }
    ch_manifest = PANGENOME_GRAPH_CONTRACT.out.manifest
    ch_variant_input = CACTUS_PANGENOME.out.all.map { taxid, products ->
        def files = products instanceof Collection ? products.toList() : [products]
        def vcf = files.find { it.name == taxid + '.gref.vcf.gz' }
        def tbi = files.find { it.name == taxid + '.gref.vcf.gz.tbi' }
        if (!vcf || !tbi) error 'Standard GREF VCF or index missing'
        tuple(taxid, vcf, tbi)
    }.join(PANGENOME_INPUT_AUDIT.out.identities).join(PANGENOME_GRAPH_CONTRACT.out.ready)
        .map { taxid, vcf, tbi, identities, check -> tuple(taxid, vcf, tbi, identities) }
    ch_versions = CACTUS_PANGENOME.out.versions
    if (params.pangenome_assignment_replay) {
        ch_assignment_input = CACTUS_PANGENOME.out.all
            .join(ch_cohorts.map { row -> tuple(row[0],row[4]) })
            .map { taxid, products, contigs ->
                def files = products instanceof Collection ? products : [products]
                def graph = files.find { it.name == taxid + '.sv.gfa.gz' }
                def paf = files.find { it.name == taxid + '.paf' }
                def gaf = files.find { it.name == taxid + '.gaf.gz' }
                if (!graph || !paf || !gaf) error 'Assignment replay requires cached Cactus SV rGFA, mapping PAF and raw GAF'
                tuple(taxid, graph, paf, gaf, contigs)
            }
    }
    ch_source_evidence = CACTUS_PANGENOME.out.all.map { taxid, products ->
        def files = products instanceof Collection ? products : [products]
        def evidence = []
        def mappings = []
        ['gaf': "${taxid}.gaf.gz", 'paf': "${taxid}.paf"].each { kind, name ->
            def selected = files.find { it.name == name }
            if (selected) {
                mappings.add(selected)
                evidence.add([kind: kind, path: selected.name])
            }
        }
        tuple(taxid, mappings, evidence)
    }
    ch_support_qc = CACTUS_PANGENOME.out.all.join(PANGENOME_GRAPH_CONTRACT.out.ready)
        .join(CACTUS_PANGENOME.out.seqfile)
        .map { taxid, products, check, seqfile ->
            def files = products instanceof Collection ? products : [products]
            tuple(taxid, files.find { it.name == taxid + '.stats.tgz' },
                  files.find { it.name == taxid + '.full.gfa.gz' },
                  file("${seqfile.parent}/versions.tsv", checkIfExists: true))
        }
    }
    emit:
    clip_stats_input = ch_clip_stats_input
    clip_sharing_input = ch_clip_sharing_input
    variant_input = ch_variant_input
    report = ch_report
    manifest = ch_manifest
    identities = PANGENOME_INPUT_AUDIT.out.identities
    versions = ch_versions
    support_qc = ch_support_qc
    source_evidence = ch_source_evidence
    assignment_input = ch_assignment_input
}
