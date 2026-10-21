nextflow.enable.dsl=2
include { CHIMERA } from '__REPO__/workflows/chimera.nf'
include { ASSEMBLY_FINISHING } from '__REPO__/workflows/assembly_finishing.nf'
include { FINALIZE_ASSEMBLY } from '__REPO__/modules/finalize_assembly.nf'
include { COLLECT_NAME_MAPS } from '__REPO__/modules/collect_name_maps.nf'
include { ASSEMBLY_RUN_SUMMARY } from '__REPO__/modules/assembly_run_summary.nf'
include { FINAL_READ_EVIDENCE } from '__REPO__/workflows/final_read_evidence.nf'
workflow {
    def root='/tmp/chimera-nxf-runtime'
    def first=[id:'sample_hap1',sample:'sample',taxid:'1',assembler:'hifiasm',hic:true]
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
    ASSEMBLY_FINISHING(CHIMERA.out.pre_finalize,
        Channel.of(tuple([sample:'sample'],file(root+'/reads.fq')),tuple([sample:'reference'],file(root+'/reads.fq'))),
        Channel.of(tuple('1','TTAGGG')), [scaffold:true,hifiasm:true])
    FINALIZE_ASSEMBLY(ASSEMBLY_FINISHING.out.assemblies)
    FINAL_READ_EVIDENCE(FINALIZE_ASSEMBLY.out.assembly,
        Channel.of(tuple([sample:'sample'],file(root+'/reads.fq')),tuple([sample:'reference'],file(root+'/reads.fq'))),
        false,Channel.empty())
    FINAL_READ_EVIDENCE.out.reference_digest.join(FINALIZE_ASSEMBLY.out.assembly).toList().view { outputs ->
        assert outputs.size()==2
        outputs.each { meta, digest, fa ->
            def actual=java.security.MessageDigest.getInstance('SHA-256').digest(fa.bytes).encodeHex().toString()
            assert digest.text.split()[0]==actual
        }
        'PASS: final HiFi mapping is bound to the finished, renamed assembly sequence'
    }
    COLLECT_NAME_MAPS(FINALIZE_ASSEMBLY.out.name_map.map { meta, nm -> nm }.collect())
    COLLECT_NAME_MAPS.out.map.view { map ->
        assert map.text.contains('sample_hap1') && map.text.contains('reference_hap1')
        if(params.fixture_manual) assert map.text.contains('s_sub_')
        'PASS: cohort name-map collection uses finalized cut identities'
    }
    ASSEMBLY_RUN_SUMMARY(Channel.value([[sample:'sample',status:'accepted',expected:1,reason:'fixture'],[sample:'reference',status:'accepted',expected:1,reason:'fixture']]),
        FINALIZE_ASSEMBLY.out.assembly.map { meta, fa -> [sample:meta.sample,id:meta.id,name:fa.name] }.collect(flat:false),
        CHIMERA.out.coordinate_report,Channel.value(file(root+'/quality.tsv')),CHIMERA.out.report)
    ASSEMBLY_RUN_SUMMARY.out.report.view { report ->
        assert report.text.contains('(assembly/chimeras/README.md)')
        assert report.text.contains('All accepted samples have their expected finalized assemblies.')
        'PASS: production run summary links the generated chimera review'
    }
    FINALIZE_ASSEMBLY.out.assembly.join(FINALIZE_ASSEMBLY.out.name_map).join(FINALIZE_ASSEMBLY.out.fai)
        .join(ASSEMBLY_FINISHING.out.teloclip_scope)
        .toList().view { outputs ->
            assert outputs.size()==2
            outputs.each { meta, fa, nm, fai, scope ->
                def lengths=fai.readLines().collectEntries { line -> def cols=line.split('\t'); [(cols[0]):cols[1].toInteger()] }
                def rows=nm.readLines().drop(1).collect { it.split('\t') }
                assert rows.size()==lengths.size()
                rows.each { row -> assert lengths[row[1]]==row[2].toInteger() }
                if(meta.id=='sample_hap1') {
                    assert lengths.size()==(params.fixture_manual ? params.fixture_cut_count+1 : 1)
                    def eligible=scope.readLines().drop(1).count { it.split('\t')[3]=='yes' }
                    assert lengths.values().sum()==(params.fixture_evidence ? 6100 : 6000)+2*eligible
                    if(params.fixture_manual) assert rows.every { it[0].startsWith('s_sub_') }
                    else assert rows[0][0]=='s'
                }
            }
            'PASS: real finishing and finalization preserve cut pieces, refreshed chromosome names and finished lengths'
        }
    CHIMERA.out.registry.view { registry ->
        assert registry.readLines().size()==3
        assert registry.text.contains('sample_hap1')
        assert registry.text.contains('reference_hap1')
        'PASS: full cohort workflow published both assemblies, including naming reference without self-PAF'
    }
    CHIMERA.out.decisions.view { decisions ->
        def lines=decisions.readLines()
        def fields=lines[0].split('\t').toList()
        def records=lines.drop(1).collect { line -> def values=line.split('\t',-1); fields.withIndex().collectEntries { name,index -> [(name):values[index]] } }
        assert records.count { it.location_kind!='alternative_region' }==((params.fixture_evidence ? 2 : 1)+(params.fixture_manual ? params.fixture_cut_count : 0))
        assert records.every { it.location_kind!='alternative_region' || (it.selected=='NO' && !it.cut_bp && it.assay_probe_bp) }
        assert decisions.text.contains('chr1') && decisions.text.contains('chr2')
        'PASS: general discovery found the synthetic chromosome change; review row emitted unselected'
    }
    CHIMERA.out.pre_finalize.toList().view { outputs ->
        assert outputs.size()==2
        if(params.fixture_manual) assert outputs.find { it[0].id=='sample_hap1' }[1].text.contains(params.fixture_cut_count==3 ? 's_sub_0_1500' : 's_sub_0_3100')
        'PASS: finishing receives both assemblies and any verified manual cut'
    }
    CHIMERA.out.report.view { report ->
        if(params.fixture_evidence) {
            def data=new groovy.json.JsonSlurper().parseText(report.parent.resolve('assemblies/sample_hap1/evidence.json').text)
            assert data.junction_policy=='junction-evidence-v1'
            assert data.events[0].junction_evidence
            assert data.events[0].junction_evidence.hifi_status=='measured'
            assert data.events[0].cut_options[0].spanning==2
            assert data.events[0].cut_options[0].spanning_strict==2
            assert data.events[0].cut_options[0].spanning_all_primary==2
            assert data.events[0].contact_profiles.size()==(params.fixture_single ? 1 : 2)
            assert data.events[0].recommendation=='REVIEW'
            assert data.events[0].localization_warning
            assert data.events[0].alternatives
            assert data.events[0].contacts.size()==(params.fixture_single ? 1 : 2)
            assert data.events[0].contacts.every { it.raw_counts.cross==0 && it.raw_counts.left_within==100 && it.raw_counts.right_within==100 && !it.calibrated }
            assert file(data.mapping.bam).exists()
        }
        if(params.fixture_manual) assert report.text.contains(params.fixture_cut_count+' actions verified')
        'PASS: source-bound read evidence, per-library measured zeros and manual application reporting'
    }
}
