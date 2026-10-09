nextflow.enable.dsl=2
include { CHIMERA } from '__REPO__/workflows/chimera.nf'
workflow {
    def root='/tmp/chimera-nxf-runtime'
    def first=[id:'sample_hap1',sample:'sample',taxid:'1',assembler:'hifiasm']
    def second=[id:'reference_hap1',sample:'reference',taxid:'1',assembler:'hifiasm']
    CHIMERA(
        Channel.of(tuple(first,file(root+'/a.fa'),file(root+'/names.tsv')),tuple(second,file(root+'/b.fa'),file(root+'/b.names.tsv'))),
        Channel.empty(),params.fixture_evidence ? Channel.of(tuple(first,file(root+'/a.agp'))) : Channel.empty(),Channel.of(tuple('1','reference_hap1')),
        Channel.of(file(root+'/quality.tsv')),
        params.fixture_evidence ? Channel.of(tuple(first,'scaffold',file(root+'/a.pairs'))) : Channel.empty(),
        params.fixture_evidence ? Channel.of(tuple(first,file(root+'/b.fa'),file(root+'/lib.tsv'))) : Channel.empty(),
        Channel.empty(),Channel.empty(),
        params.fixture_evidence ? Channel.of(tuple([sample:'sample'],file(root+'/reads.fq'))) : Channel.empty(),
        [harmonize:true,scaffold:false])
    CHIMERA.out.registry.view { registry ->
        assert registry.readLines().size()==3
        assert registry.text.contains('sample_hap1')
        assert registry.text.contains('reference_hap1')
        'PASS: full cohort workflow published both assemblies, including naming reference without self-PAF'
    }
    CHIMERA.out.decisions.view { decisions ->
        assert decisions.readLines().size()==(params.fixture_manual ? 3 : 2)
        assert decisions.text.contains('chr1') && decisions.text.contains('chr2')
        'PASS: general discovery found the synthetic chromosome change; review row emitted unselected'
    }
    CHIMERA.out.pre_finalize.toList().view { outputs ->
        assert outputs.size()==2
        if(params.fixture_manual) assert outputs.find { it[0].id=='sample_hap1' }[1].text.contains('s_sub_0_3100')
        'PASS: finishing receives both assemblies and any verified manual cut'
    }
    CHIMERA.out.report.view { report ->
        if(params.fixture_evidence) {
            def data=new groovy.json.JsonSlurper().parseText(report.parent.resolve('assemblies/sample_hap1/evidence.json').text)
            assert data.events[0].cut_options[0].spanning==2
            assert data.events[0].recommendation=='RETAIN'
            assert data.events[0].contacts.size()==2
            assert data.events[0].contacts.every { it.raw_counts.cross==0 && it.raw_counts.left_within==100 && it.raw_counts.right_within==100 && !it.calibrated }
            assert file(data.mapping.bam).exists()
        }
        if(params.fixture_manual) assert report.text.contains('1 actions verified')
        'PASS: source-bound read evidence, per-library measured zeros and manual application reporting'
    }
}
