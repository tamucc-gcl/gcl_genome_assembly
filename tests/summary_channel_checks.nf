nextflow.enable.dsl=2
include { graphReportPath; rebaseReportMarkdown } from '__REPO__/functions/report_artifacts.nf'
workflow {
    def fixtureRepo = '__REPO__'
    def root=file('/tmp/chimera-nxf-runtime/summary-fragments');root.mkdirs()
    def firstFragment=file("${root}/pangenome_build_report.md");firstFragment.text='FIRST_FRAGMENT'
    def secondFragment=file("${root}/pangenome_input_report.md");secondFragment.text='SECOND_FRAGMENT'
    def mode=System.getenv('SUMMARY_CASE') ?: 'multiple'
    ch_reports=mode=='empty' ? Channel.empty() : Channel.of(['2',secondFragment],['1',firstFragment])
__PANGENOME_REPORT_COLLECTION__
    ch_pangenome_report_for_report.toList().view { paths ->
        assert paths.size()==1
        if(mode=='empty') {
            assert paths[0].name=='NO_PANGENOME'
            assert paths[0].exists()
        } else {
            def body=paths[0].text
            assert body.indexOf('FIRST_FRAGMENT')>=0
            assert body.indexOf('SECOND_FRAGMENT')>body.indexOf('FIRST_FRAGMENT')
        }
        'PASS: single report input for '+mode
    }
}
