nextflow.enable.dsl=2
include { CHECK } from '__REPO__/scripts/comparisons/fcs_adaptor_acceptance.nf'
workflow {
    def root=file('/tmp/fcs-check-staging');root.mkdirs()
    def original=[]
    def cleaned=[]
    def actions=['#accession\tlength\taction\trange\tname']
    for(i in 1..10) {
        def id=String.format('seq_%05d',i)
        def sequence=i<=5 ? 'A'*1000 : 'ACGT'*25100
        original.add('>'+id+'\n'+sequence+'\n')
        if(i<=5) {
            actions.add("${id}\t${sequence.size()}\tACTION_TRIM\t1..58\tfixture")
            cleaned.add('>'+id+'~59..'+sequence.size()+'\n'+sequence.substring(58)+'\n')
        } else {
            actions.add("${id}\t${sequence.size()}\tACTION_TRIM\t100001..100058\tfixture")
            cleaned.add('>'+id+'~1..100000\n'+sequence.substring(0,100000)+'\n')
            cleaned.add('>'+id+'~100059..'+sequence.size()+'\n'+sequence.substring(100058)+'\n')
        }
    }
    original.add('>seq_00018\n'+('A'*522)+'\n')
    actions.add('seq_00018\t522\tACTION_EXCLUDE\t\tfixture')
    def source=root.resolve('source');source.mkdirs()
    def clean=root.resolve('cleaned');clean.mkdirs()
    // Deliberately identical basenames for both FASTAs and both reports.
    source.resolve('assembly.fasta').text=original.join('')
    clean.resolve('assembly.fasta').text=cleaned.join('')
    def initial=root.resolve('initial');initial.mkdirs()
    def repeat=root.resolve('repeat');repeat.mkdirs()
    initial.resolve('ncbi_fixture.fcs_adaptor_report.txt').text=actions.join('\n')+'\n'
    repeat.resolve('ncbi_fixture.fcs_adaptor_report.txt').text=actions[0]+'\n'
    CHECK(source.resolve('assembly.fasta'),initial.resolve('ncbi_fixture.fcs_adaptor_report.txt'),
        clean.resolve('assembly.fasta'),repeat.resolve('ncbi_fixture.fcs_adaptor_report.txt'),
        ['tests/fcs_adaptor_acceptance.py','py_scripts/contamination_review.py','py_scripts/contamination_cohort_report.py']
            .collect { file('__REPO__/'+it,checkIfExists:true) })
    CHECK.out.view { result ->
        def data=new groovy.json.JsonSlurper().parse(result)
        assert data.status=='PASS'
        assert data.terminal_trims==5 && data.internal_splits==5 && data.excluded_sequences==1
        assert data.removed_bp==1102 && data.repeat_adaptor_actions==0
        'PASS real checker with identical input basenames and distinct initial/repeat reports'
    }
}
