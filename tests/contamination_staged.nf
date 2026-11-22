nextflow.enable.dsl = 2

// Exercise real search/report channel wiring. A tiny database and tools are
// supplied by the caller; this test does not download production databases.
params.fixture_fasta = null
params.fixture_empty = null
params.fixture_taxonomy = null
params.fixture_protein = null
params.fixture_nt = ''
params.fixture_mode = 'staged'
include { CONTAMINATION_SEARCH } from '../workflows/contamination_search.nf'

workflow {
    def sources = Channel.of(tuple('populated', file(params.fixture_fasta), 7955),
                            tuple('empty', file(params.fixture_empty), 7955))
    CONTAMINATION_SEARCH(sources,
        Channel.of(tuple('populated', [], [], []), tuple('empty', [], [], [])),
        params.fixture_taxonomy,
        [protein: params.fixture_protein, stage: 'post_fcs', window: 100000,
         chunk: 5000000, mode: params.fixture_mode, nt: params.fixture_nt,
         mega: !!params.fixture_nt, targeted: false, followup: 100000, rank: 'phylum', nested: true],
        file("${projectDir}/.."), [])
    CONTAMINATION_SEARCH.out.review.toList().view { rows ->
        assert rows.size() == 2
        rows.each { id, review ->
            def evidence = new groovy.json.JsonSlurper().parseText(review.resolve('evidence.json').text)
            assert evidence.search_mode == params.fixture_mode
            assert evidence.protein_searched_bp <= evidence.source_bp
            assert evidence.megablast_searched_bp <= evidence.source_bp
            if (id == 'empty') assert evidence.search_queries == 0 && evidence.protein_searched_bp == 0
            else assert evidence.search_queries > 0 && evidence.protein_searched_bp > 0
            if (params.fixture_mode == 'staged') assert review.resolve('sampling_windows.tsv').size() > 0
        }
        "${params.fixture_mode}: populated and zero-query assemblies reached reports"
    }
}
