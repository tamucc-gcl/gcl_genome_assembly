nextflow.enable.dsl = 2

include { MAPPING_QC as POST_FCS_MAPPING } from '../modules/mapping_qc.nf'
include { CONTAMINATION_SEARCH } from './contamination_search.nf'
include { CONTAMINATION_DECISIONS; APPLY_CONTAMINATION_REVIEW } from '../modules/contamination_cleanup.nf'
include { DOWNLOAD_CONTAMINATION_NT; VERIFY_CONTAMINATION_NT } from '../modules/download_contamination_nt.nf'
include { CONTAMINATION_PROTEIN_DB } from './contamination_protein_db.nf'
include { DOWNLOAD_CONTAMINATION_PROTEIN_TAXONOMY as REUSE_CONTAMINATION_PROTEIN_TAXONOMY } from '../modules/contamination_protein_source.nf'

workflow CONTAMINATION_POST_FCS {
    take:
    assemblies
    reads
    taxdump

    main:
    def base = file("${moduleDir}/..").normalize()
    def mapping = Channel.empty()
    def mappingVersions = Channel.empty()
    def coverage = assemblies.map { meta, fasta -> tuple(meta.id, [], [], []) }
    def reviews = Channel.empty()
    if (params.run_contamination_evidence) {
        def taxonomy = taxdump.collect().map { dirs -> dirs[0].toString() }
        def protein = params.contamination_protein_db_dir
        if (params.contamination_protein_db_root) {
            def source = params.contamination_protein_source
            if (source && !source.toString().startsWith('https://'))
                source = file(source).toAbsolutePath().toString()
            CONTAMINATION_PROTEIN_DB(
                [root: file(params.contamination_protein_db_root).toAbsolutePath().toString(),
                 release: params.contamination_protein_release.toString(), source: source,
                 sha256: params.contamination_protein_sha256, format: params.contamination_protein_format,
                 taxonomy: params.contamination_protein_taxonomy_dir ? file(params.contamination_protein_taxonomy_dir).toAbsolutePath().toString() : null],
                base)
            protein = CONTAMINATION_PROTEIN_DB.out.directory
            taxonomy = CONTAMINATION_PROTEIN_DB.out.taxdump
        } else if (params.contamination_protein_taxonomy_dir) {
            REUSE_CONTAMINATION_PROTEIN_TAXONOMY(file(params.contamination_protein_taxonomy_dir).toAbsolutePath().toString())
            taxonomy = REUSE_CONTAMINATION_PROTEIN_TAXONOMY.out.directory
        }
        def assemblyInput = assemblies.map { meta, fasta -> tuple(meta.id, fasta, meta.taxid as int) }
        if (params.contamination_read_support) {
            def mappingInput = assemblies.map { meta, fasta -> tuple(meta.sample, meta, fasta) }
                .combine(reads.map { meta, reads_for_qc -> tuple(meta.sample, [reads_for_qc]) }, by: 0)
                .map { sample, meta, fasta, wrapped -> tuple(meta, fasta, wrapped[0]) }
            POST_FCS_MAPPING(mappingInput, true)
            mapping = POST_FCS_MAPPING.out.results
            mappingVersions = POST_FCS_MAPPING.out.versions
            coverage = POST_FCS_MAPPING.out.alignment.map { meta, bam, bai -> tuple(meta.id, bam, bai, []) }
        }
        def nt = params.contamination_nucleotide_db ? file(params.contamination_nucleotide_db).toAbsolutePath().toString() : ''
        def ntManifest = []
        if (params.contamination_nt_db_root) {
            def root = file(params.contamination_nt_db_root).toAbsolutePath().toString()
            def downloader = file(base.resolve('py_scripts/download_contamination_nt.py'), checkIfExists: true)
            DOWNLOAD_CONTAMINATION_NT(root, 'nt', params.contamination_nt_release, downloader)
            VERIFY_CONTAMINATION_NT(DOWNLOAD_CONTAMINATION_NT.out.bundle, root,
                'nt', params.contamination_nt_release, downloader)
            nt = VERIFY_CONTAMINATION_NT.out.prefix.collect().map { prefixes -> prefixes[0] }
            ntManifest = VERIFY_CONTAMINATION_NT.out.manifest
        }
        def mega = !!(params.contamination_nucleotide_db || params.contamination_nt_db_root)
        CONTAMINATION_SEARCH(assemblyInput, coverage, taxonomy,
            [protein: protein, stage: 'post_fcs', window: 100000, mode: params.contamination_search_mode,
             chunk: 5000000, nt: nt, mega: mega, targeted: params.contamination_targeted_blastn,
             followup: 100000, rank: params.contamination_host_rank, nested: true],
            base, ntManifest)
        reviews = CONTAMINATION_SEARCH.out.review
    }

    def cleaned = assemblies
    if (params.contamination_review) {
        def sheet = file(params.contamination_review, checkIfExists: true)
        def helper = file(base.resolve('py_scripts/contamination_review.py'), checkIfExists: true)
        CONTAMINATION_DECISIONS(sheet, assemblies.map { meta, fasta -> meta.id }.toList(), helper)
        def decisionPlan = CONTAMINATION_DECISIONS.out.collect().map { files ->
            [selected: new groovy.json.JsonSlurper().parseText(files[0].text)]
        }
        def routed = assemblies.combine(decisionPlan).branch { meta, fasta, plan ->
            apply: meta.id.toString() in plan.selected
                return tuple(meta, fasta)
            keep: true
                return tuple(meta, fasta)
        }
        APPLY_CONTAMINATION_REVIEW(routed.apply, sheet, helper)
        cleaned = routed.keep.mix(APPLY_CONTAMINATION_REVIEW.out.assembly)
    }

    emit:
    assembly = cleaned
    review = reviews
    mapping_qc = mapping
    versions = mappingVersions
}
