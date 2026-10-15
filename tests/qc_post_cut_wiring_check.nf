nextflow.enable.dsl=2
include { ASSEMBLY_QC } from '__REPO__/workflows/assembly_qc.nf'

workflow CHIMERA {
    main:
    def root = '/tmp/chimera-nxf-runtime'
    def meta = { sample, hap, count -> [id:"${sample}_${hap}",sample:sample,haplotype:hap,n_hap:count,taxid:'1',hifi:true] }
    all = Channel.of(
        tuple(meta('changed','hap2',2),file(root+'/b.fa'),file(root+'/names.tsv')),
        tuple(meta('unchanged','hap1',2),file(root+'/a.fa'),file(root+'/names.tsv')),
        tuple(meta('changed','hap1',2),file(root+'/a.fa'),file(root+'/names.tsv')),
        tuple(meta('unchanged','hap2',2),file(root+'/a.fa'),file(root+'/names.tsv')),
        tuple(meta('singleton','primary',1),file(root+'/a.fa'),file(root+'/names.tsv')))
    cut = all.filter { m,f,n -> (m.sample == 'changed' && m.haplotype == 'hap1') || m.sample == 'singleton' }
    emit:
    broken = cut
    pre_finalize = all
}

workflow {
    def root = '/tmp/chimera-nxf-runtime'
    CHIMERA()
    def run_all_qc = true
    def chimera_on = true
    ch_staged_assemblies = Channel.empty()
__POST_CUT_QC__
    ch_staged_assemblies.toList().view { rows ->
        assert rows.size() == 3
        assert rows.collect { it[0].id }.sort() == ['changed_hap1','changed_hap2','singleton_primary']
        assert rows.every { it[1] == 'chimera_broken' }
        'PASS: post-cut checkpoint includes uncut sister, excludes unchanged sample, supports singleton'
    }
    ASSEMBLY_QC(ch_staged_assemblies.map { m,s,f -> tuple(m,f) },
        Channel.of(tuple([sample:'changed'],file(root+'/reads.fq')),tuple([sample:'singleton'],file(root+'/reads.fq'))),
        Channel.of(tuple([sample:'changed'],file(root+'/names.tsv')),tuple([sample:'singleton'],file(root+'/names.tsv'))),
        Channel.value(['1':'stub']), 'chimera_broken')
    ASSEMBLY_QC.out.assembly_summary.toList().view { rows ->
        assert rows.size() == 2
        assert rows.collect { it[0] }.sort() == ['changed','singleton']
        'PASS: production QC grouping completes paired and singleton summaries'
    }
}
