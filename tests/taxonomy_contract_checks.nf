nextflow.enable.dsl=2
include { RESOLVE_TAXONOMY } from '__REPO__/modules/resolve_taxonomy.nf'
workflow {
    def root=file('/tmp/chimera-nxf-runtime/taxonomy-fixtures');root.mkdirs()
    def db=file("${root}/taxdump");db.mkdirs()
    file("${db}/names.dmp").text='fixture\n'
    file("${db}/nodes.dmp").text='fixture\n'
    RESOLVE_TAXONOMY(Channel.of(tuple('1',db)))
    RESOLVE_TAXONOMY.out.tsv.view { taxid,tsv ->
        def row=tsv.text.trim().readLines().last().split('\t')
        assert row.size()==10 && row[0]=='1' && row[9].contains('Metazoa')
        'PASS: taxonomy resolution emits a complete species-scoped record'
    }
}
