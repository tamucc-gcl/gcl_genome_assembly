nextflow.enable.dsl = 2

params.fixture_fasta = null
params.fixture_taxonomy = null
include { CONTAMINATION_POST_FCS } from '../workflows/contamination_post_fcs.nf'

workflow {
    def label = 'empty'
    CONTAMINATION_POST_FCS(
        Channel.of(tuple([id: "fixture_${label}_hap1", sample: 'fixture', taxid: 7955], file(params.fixture_fasta))),
        Channel.empty(), Channel.value(file(params.fixture_taxonomy)))
    CONTAMINATION_POST_FCS.out.review.toList().view { rows ->
        assert rows.size() == 1
        def evidence = new groovy.json.JsonSlurper().parseText(rows[0][1].resolve('evidence.json').text)
        assert evidence.source_bp == 200 && evidence.search_queries == 0 && evidence.chunks == 0
        assert evidence.protein_database.protein_records == 2
        'Managed database reached the real post-FCS evidence/report workflow; zero-query assembly preserved'
    }
    CONTAMINATION_POST_FCS.out.assembly.toList().view { rows ->
        assert rows.size() == 1 && rows[0][1].toRealPath() == file(params.fixture_fasta).toRealPath()
        'Unreviewed input preserved'
    }
}
