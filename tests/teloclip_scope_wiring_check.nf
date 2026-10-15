nextflow.enable.dsl=2
include { ASSEMBLY_FINISHING } from '__REPO__/workflows/assembly_finishing.nf'
workflow {
    def root='/tmp/chimera-nxf-runtime'
    def inputs=['mixed','fragmented','unclassified'].collect { sample ->
        tuple([id:sample,sample:sample,taxid:'1',assembler:'hifiasm',hic:false],
              file(root+"/${sample}.fa"),
              file(root+(sample=='unclassified' ? '/NO_HARMONIZE' : "/${sample}.tsv"))) }
    ASSEMBLY_FINISHING(Channel.fromList(inputs),
        Channel.fromList(inputs.collect { row -> tuple([sample:row[0].sample],file(root+'/reads.fq')) }),
        Channel.of(tuple('1','TTAGGG')), [scaffold:false,hifiasm:true])
    ASSEMBLY_FINISHING.out.assemblies.toList().view { rows ->
        assert rows.size()==3
        rows.each { m,fa,nm ->
            def input=file(root+"/${m.id}.fa").text
            if(m.id=='mixed') {
                assert fa.text==input.replace('AAAACCCC','AAAACCCCTT')
                assert fa.text.contains('>small\nGGGG\n>uncertain\nTTTT\n')
            } else {
                assert fa.text==input
                assert !fa.parent.resolve("${m.id}.teloclip_overhangs.bam").exists()
                assert fa.parent.resolve('versions.tsv').text.contains('skipped_no_chromosome_scale_scaffolds')
            }
            def audit=fa.parent.resolve("${m.id}.teloclip_scope.tsv").readLines()
            assert audit.size()==4
            assert audit.drop(1).count { it.split('\t')[3]=='yes' }==(m.id=='mixed' ? 1 : 0)
        }
        'PASS: production finishing extends only the classified scaffold; fragmented/unclassified assemblies skip mapping and preserve every sequence'
    }
}
