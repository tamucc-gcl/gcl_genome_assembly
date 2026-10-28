nextflow.enable.dsl=2
include { RUN_PROVENANCE } from '__REPO__/workflows/run_provenance.nf'
workflow {
    def root=file('/tmp/chimera-nxf-runtime/provenance-fixtures');root.mkdirs()
    def one=file("${root}/one/versions.tsv");one.parent.mkdirs()
    def two=file("${root}/two/versions.tsv");two.parent.mkdirs()
    one.text='tool\tversion\nminimap2\t2.29\n'
    two.text='process\ttool\tversion\nMAP\tminimap2\t2.30\n'
    def mode=System.getenv('PROVENANCE_CASE') ?: 'nonempty'
    RUN_PROVENANCE(mode=='empty' ? Channel.empty() : Channel.of(one,two))
    RUN_PROVENANCE.out.versions.toList().view { paths ->
        assert paths.size()==1
        assert paths[0].text.readLines().size()==(mode=='empty' ? 1 : 3)
        assert paths[0].text.startsWith('tool\tversion\tprocess\tsource')
        'PASS: always-on provenance closes for '+mode+' versions with real production collector'
    }
}
