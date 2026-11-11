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
    graph=mode=='multispecies' ? Channel.fromList(['sharing_summary.tsv','sharing_report.md','source_sequences.tsv','source_audit.json','source_report.md','pangenome_manifest.tsv'].collect { tuple('1',file("${root}/${it}")) }) : Channel.empty()
    snails=mode!='disabled' ? Channel.fromList(data.assemblies.collect { tuple(it.id,'final',file("${root}/${it.id}_final_snail.svg")) }) : Channel.empty()
    annotations=mode!='disabled' ? Channel.value(tuple([sample:'sample'],file("${root}/annotation.gb"))) : Channel.empty()
    psmc=mode=='psmc' ? Channel.fromList(['README.md','summary.tsv','combined_history.png'].collect { file("${root}/psmc/${it}") }) : Channel.empty()
    def reviewMode=mode in ['uncut_review','single_cut_review','multi_cut_review']
    def decisionFile=reviewMode ? file("${root}/review.tsv") : no
    def verificationFiles=mode=='single_cut_review' ? [file("${root}/sample_hap1.verification.json")] :
        mode=='multi_cut_review' ? [file("${root}/sample_hap1.verification.json"),file("${root}/sample_hap2.verification.json")] :
        [file("${root}/NO_PAIRS")]
    MAIN_REPORTING(eligibility,expected,names,REPORT_PAINTINGS.out.blocks,qc,
        Channel.value(decisionFile),Channel.value(verificationFiles),graph,Channel.empty(),
        Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),
        annotations,Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),snails,
        Channel.value(data.assemblies.collect { [sample:it.sample,hifi:false,hic:false,shortread:false] }.unique { it.sample }),Channel.empty(),
        Channel.value([file("${root}/software_versions.tsv"),file("${root}/run_provenance.json")]),psmc)
    MAIN_REPORTING.out.report.view { report ->
        assert report.text.contains('| sample | Finalized |')
        assert report.text.contains('## Chromosome structure')==(mode!='disabled')
        assert report.text.contains('## Assembly quality')==(mode!='disabled')
        assert report.text.contains('sample_hap1_final_snail.svg')==(mode!='disabled')
        assert !file("${report.parent}/reports/main/report_completeness.tsv").text.contains('\tMISSING\t')
        assert file("${report.parent}/reports/main/artifact_registry.tsv").text.contains('provenance/software_versions.tsv')
        assert report.text.contains('## Pangenome')==(mode=='multispecies')
        assert report.text.contains('[Download the pangenome graph (CLIP GBZ)](pangenome/1/1.gbz)')==(mode=='multispecies')
        assert !report.text.contains('pangenome/1/1.gref.gbz')
        assert report.text.contains('## Demographic history (PSMC)')==(mode=='psmc')
        assert report.text.contains('![PSMC demographic histories for fitted individuals](post_assembly/psmc/combined_history.png)')==(mode=='psmc')
        if(mode=='psmc') assert report.text.contains('| sample | assemblies | completed | 7.50 | 75.00% | 1 |')
        if(reviewMode) {
            assert report.text.contains('## Assembly decisions and changes')
            assert report.text.contains(mode=='uncut_review' ? 'Review pending' : 'Cut applied and verified')
            def spec=new groovy.json.JsonSlurper().parse(file("${report.parent}/report_spec.json"))
            assert spec.verifications.size()==(mode=='multi_cut_review' ? 2 : 1)
            assert spec.verifications.every { it.startsWith('verification/') && file("${report.parent}/${it}").isFile() }
        }
        if(mode=='multispecies') {
            assert file("${report.parent}/reports/main/needs_attention.tsv").text.contains('sample_hap1 / chr1_1')
            assert !file("${report.parent}/reports/main/needs_attention.tsv").text.contains('sample_hap2 / chr1_1')
            assert report.text.contains('Distinct graph sequence (Mb)')
            assert report.text.contains('All 2 individuals')
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
