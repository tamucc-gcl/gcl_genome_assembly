nextflow.enable.dsl=2
include { FINAL_VIZ } from '__REPO__/workflows/final_viz.nf'
include { PAIRWISE_COMPARISONS } from '__REPO__/workflows/pairwise_comparisons.nf'
workflow {
    def mode=System.getenv('PAIRWISE_CASE') ?: 'within_sample'
    def root=file('/tmp/chimera-nxf-runtime/comparison-fixtures');root.mkdirs()
    def records=[['first','alpha','1'],['second','alpha','1'],['single','beta','1'],
                 ['third','gamma','2'],['fourth','gamma','2'],['lonely','delta','2']].collect { id,sample,taxid ->
        def fa=file("${root}/${id}.fasta")
        def fai=file("${root}/${id}.fasta.fai")
        if (System.getenv('PAIRWISE_RESUME_CHECK')!='true') {
            fa.text='>chr1\n'+'A'*10000+'\n'
            fai.text='chr1\t10000\t6\t10000\t10001\n'
        }
        tuple([id:id,sample:sample,taxid:taxid],fa,fai)
    }
    if (mode=='no_pairs') records=records.findAll { it[0].id in ['single','lonely'] }
    if (mode=='duplicate_id') records[1][0].id='first'
    if (mode=='conflicting_species') records[1][0].taxid='2'
    if (mode=='wrong_dictionary') records[0][2].text='chr1\t9999\t6\t9999\t10000\n'
    def expected=mode in ['empty','no_pairs','disabled'] ? 0 : (mode=='all_species' ? 6 : 2)
    if (mode=='empty') {
        PAIRWISE_COMPARISONS(Channel.empty())
        PAIRWISE_COMPARISONS.out.summary.toList().view { paths ->
            assert paths.size()==1 && paths[0].text.readLines().size()==1
            'PASS: empty cohort produces an explicit header-only comparison summary'
        }
    } else {
        FINAL_VIZ(Channel.fromList(records),Channel.fromList(records.collect { m,fa,fai -> tuple(m.id,fa,'CCCTAA') }))
        FINAL_VIZ.out.dotplot.toList().view { rows ->
            assert rows.size()==expected
            assert rows.every { it[2].size()>1000 }
            'PASS: correct number of actual rendered dotplots'
        }
        FINAL_VIZ.out.riparian.toList().view { rows ->
            assert rows.size()==expected && rows.every { it[2].size()>1000 }
            'PASS: correct number of actual rendered ribbon plots'
        }
        FINAL_VIZ.out.pairwise_summary.toList().view { paths ->
            assert paths.size()==1
            if (mode=='disabled') assert paths[0].name=='NO_PAIRWISE'
            else {
                def lines=paths[0].text.readLines();assert lines.size()==expected+1
                def pairs=lines.drop(1).collect { it.split('\t')[0..1] as Set }
                assert !pairs.any { pair -> pair.any { it in ['first','second','single'] } && pair.any { it in ['third','fourth','lonely'] } }
                if (mode=='within_sample') assert pairs as Set==[['first','second'] as Set,['third','fourth'] as Set] as Set
                if (mode=='no_alignments') assert lines.drop(1).every { it.split('\t')[5]=='0' && it.split('\t')[8]=='NA' }
            }
            'PASS: species-aware pairing and final-coordinate QC summary'
        }
        FINAL_VIZ.out.telomere_summary.toList().view { paths ->
            assert paths.size()==1 && paths[0].text.readLines().size()==records.size()+1
            'PASS: complete FINAL_VIZ telomere/QUAST routing'
        }
    }
}
