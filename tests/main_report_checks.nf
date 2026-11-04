nextflow.enable.dsl=2
include { REPORT_PAINTINGS; MAIN_REPORTING } from '__REPO__/workflows/main_report.nf'
include { graphReportPath; rebaseReportMarkdown } from '__REPO__/functions/report_artifacts.nf'
workflow {
    def root=file('/tmp/chimera-nxf-runtime/main-report-fixtures')
    def mode=System.getenv('MAIN_REPORT_CASE') ?: 'enabled'
    eligibility=Channel.value(file("${root}/eligibility.json"))
    REPORT_PAINTINGS(eligibility)
    def data=new groovy.json.JsonSlurper().parse(file("${root}/eligibility.json"))
    expected=Channel.value(data.assemblies.groupBy { it.sample }.collect { sample,rows ->
        [sample:sample,status:'accepted',expected:rows.size(),reason:'Accepted input'] })
    names=Channel.fromList(data.assemblies.collect { row -> [[id:row.id],file("${root}/${row.id}.name_map.tsv")] })
    def no=file('__REPO__/assets/NO_REPORT_INPUT')
    qc=Channel.empty()
    if(mode!='disabled') qc=Channel.value(file("${root}/qc.csv"))
    MAIN_REPORTING(eligibility,expected,names,REPORT_PAINTINGS.out.blocks,qc,
        Channel.value(no),Channel.value([no]),Channel.empty(),Channel.empty(),
        Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),
        Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),Channel.value(no))
    MAIN_REPORTING.out.report.view { report ->
        assert report.text.contains('| sample | Finalized |')
        assert report.text.contains('## Chromosome structure')==(mode!='disabled')
        assert report.text.contains('## Assembly quality')==(mode!='disabled')
        assert !report.text.contains('## Pangenome')
        if(mode=='multispecies') {
            assert report.text.contains('chromosome_structure_1.png')
            assert report.text.contains('chromosome_structure_2.png')
            def scope=file("${report.parent}/reports/main/display_scope.tsv").readLines().drop(1)
            assert scope.size()==4
            assert scope.findAll { it.startsWith('sample_B_') }.every { it.endsWith('\t2') }
        }
        'PASS: production main report '+mode+' without post-assembly analyses'
    }
    assert graphReportPath('1','sharing_figures.md')=='pangenome/1/figures/sharing_figures.md'
    def rebased=rebaseReportMarkdown('# Title\n![Image](plots/a.png)\n[Data](../data.tsv)\n<img src="x.png" />','pangenome/1/representation/report.md')
    assert rebased.contains('#### Title')
    assert rebased.contains('(pangenome/1/representation/plots/a.png)')
    assert rebased.contains('(pangenome/1/data.tsv)')
    assert rebased.contains('src="pangenome/1/representation/x.png"')
}
