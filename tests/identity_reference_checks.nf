nextflow.enable.dsl=2
include { INPUT_TAXONOMY } from '__REPO__/workflows/input_taxonomy.nf'
include { FIND_MITO_REFERENCE } from '__REPO__/modules/find_mito_reference.nf'
workflow {
    def mode=System.getenv('ORG_CASE') ?: 'identity_mixed'
    if (mode.startsWith('reference_')) {
        FIND_MITO_REFERENCE(Channel.of(tuple('11','Species one')))
        FIND_MITO_REFERENCE.out.ref_info.toList().view { rows ->
            assert rows.size()==1 && rows[0][1].text.contains('ACC\tSpecies one\t4\tSpecies one')
            'PASS: validated mitochondrial reference pair'
        }
    } else {
        def root=file('/tmp/chimera-nxf-runtime/identity-fixtures');root.mkdirs()
        def dump=file("${root}/db");dump.mkdirs()
        file("${dump}/names.dmp").text='11\t|\tSpecies one\t|\t\t|\tscientific name\t|\n22\t|\tSpecies two\t|\t\t|\tscientific name\t|\n33\t|\tSpecies two\t|\t\t|\tcommon name\t|\n'+
            (mode=='identity_ambiguous' ? '33\t|\tSpecies two\t|\t\t|\tscientific name\t|\n' : '')
        file("${dump}/nodes.dmp").text='11\t|\t1\t|\tspecies\t|\n'+
            (mode=='identity_missing_node' ? '' : '22\t|\t1\t|\tspecies\t|\n33\t|\t1\t|\tspecies\t|\n')
        def known=[sample:'known',taxid:'11',species:'Explicit name override',assembler:'hifiasm']
        def named=[sample:'named',taxid:null,species:(mode=='identity_unknown' ? 'Unknown species' : 'species   TWO'),assembler:'spades']
        if (mode=='identity_missing') named.species=null
        if (mode=='identity_invalid_taxid') named.taxid='11;bad'
        def records=[tuple(known,[file:'known.fq'])]
        if (mode!='identity_known_only') records << tuple(named,[file:'named.fq'])
        INPUT_TAXONOMY(Channel.fromList(records),Channel.value(dump))
        INPUT_TAXONOMY.out.samples.toList().view { rows ->
            assert rows.size()==(mode=='identity_known_only' ? 1 : 2)
            assert rows.find { it[0].sample=='known' }==tuple(known,[file:'known.fq'])
            if (mode!='identity_known_only') {
                assert rows.find { it[0].sample=='named' }[0].taxid=='22'
                assert rows.find { it[0].sample=='named' }[0].species=='Species two'
                assert rows.find { it[0].sample=='named' }[1]==[file:'named.fq']
            }
            'PASS: identity resolution preserves every sample and numeric metadata'
        }
    }
}
