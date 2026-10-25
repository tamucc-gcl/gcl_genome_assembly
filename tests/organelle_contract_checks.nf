nextflow.enable.dsl=2
include { ORGANELLE } from '__REPO__/workflows/organelle.nf'

workflow {
    def mode = System.getenv('ORG_CASE') ?: 'circular'
    def successCase=mode in ['circular','linear','incomplete','no_recovery','retry','fungus_retry','multi_target','multi_species']
    def hifiCase=mode in ['missing_reference','hifi_missing_gcode']
    def root=file('/tmp/chimera-nxf-runtime/organelle-fixtures'); root.mkdirs()
    def r1=file("${root}/r1.fq");r1.text='@read\nACGT\n+\nIIII\n'
    def r2=file("${root}/r2.fq");r2.text=r1.text
    def org=(mode=='fungus_retry') ? 'fungus_mt' : 'animal_mt'
    def specs=(mode=='multi_target') ? ['embplant_pt','embplant_mt'] : [org]
    def targets=specs.collect { [type:it,coverage:50,recursion:10,kmers:'21,45'] }
    def reads=['alpha','beta'].collect { name ->
        tuple([sample:name,taxid:(mode=='multi_species' && name=='beta' ? '2' : '1'),hifi:hifiCase,mito_tool:(hifiCase ? 'mitohifi' : 'mitofinder')],r1,r1,r2)
    }
    def targetContext=(mode=='missing_targets') ? Channel.empty() : Channel.of(tuple('1',targets))
    if (mode=='duplicate_targets') targetContext=Channel.of(tuple('1',targets),tuple('1',targets))
    if (mode=='multi_species') targetContext=Channel.of(tuple('1',targets),tuple('2',[[type:'fungus_mt',coverage:50,recursion:10,kmers:'21,45']]))
    def codes=(mode in ['missing_gcode','hifi_missing_gcode']) ? Channel.empty() : Channel.of(tuple('1',2))
    if (mode=='multi_species') codes=Channel.of(tuple('1',2),tuple('2',4))
    def references=(mode=='hifi_missing_gcode') ? Channel.of(tuple('1',r1,r2)) : Channel.empty()
    ORGANELLE(Channel.fromList(reads),references,codes,targetContext,
        [hifi:hifiCase,shortread_organelle:!hifiCase])
    ORGANELLE.out.baits.toList().view { rows ->
        if (successCase) {
            assert rows*.getAt(0).sort()==['alpha','beta']
            assert rows.every { it[1].size()==(mode=='no_recovery' ? 0 : specs.size()) }
        }
        'PASS: each short-read sample releases exactly one correctly sized bait set'
    }
    ORGANELLE.out.annotation.toList().view { rows ->
        if (successCase) {
            assert rows.size()==(mode in ['no_recovery','multi_target'] ? 0 : 2)
            assert rows.every { it[1].text.contains('LOCUS') && it[1].text.contains('/gene="cox1"') }
        }
        'PASS: actual GenBank converter and annotation output contract'
    }
    ORGANELLE.out.notes.toList().view { rows ->
        if (successCase) assert rows.size()==(mode=='multi_target' ? 4 : 0)
        'PASS: plant organelles retain explicit annotation notes'
    }
    ORGANELLE.out.recovery_status.toList().view { rows ->
        if (successCase) assert rows.size()==2*specs.size()
        def expected=mode=='no_recovery' ? 'no_recovery' : (mode in ['linear','incomplete'] ? 'linear' : 'circular')
        if (successCase) assert rows.every { it[2].text.trim().readLines().last().split('\t')[2]==expected }
        if (mode=='multi_species') assert rows.every { it[1]==(it[0].taxid=='1' ? 'animal_mt' : 'fungus_mt') }
        'PASS: recovery states distinguish circular, linear and completed no-recovery'
    }
    ORGANELLE.out.stats.toList().view { rows ->
        if (mode=='multi_species') {
            assert rows.size()==2
            assert rows.every { it[1].text.trim().readLines().last().split('\t')[6]==(it[0].taxid=='1' ? '2' : '4') }
        }
        'PASS: annotation genetic codes stay attached to their own species'
    }
}
