nextflow.enable.dsl=2
include { DECONTAMINATE_ASSEMBLY as CLEAN } from '__REPO__/workflows/decontaminate_assembly.nf'
include { FCS_FINAL_SCREEN as AUDIT } from '__REPO__/workflows/fcs_final_screen.nf'
include { FCS_DIAGNOSTICS as PLOTS } from '__REPO__/workflows/fcs_diagnostics.nf'
include { CONTAMINATION_REPORTING as REPORT } from '__REPO__/workflows/contamination_reporting.nf'
include { FCS_OVERVIEW as OVERVIEW } from '__REPO__/workflows/fcs_overview.nf'
params.outdir='/tmp/fcs-workflow-check/results'
workflow {
    def root=file('/tmp/fcs-workflow-check/inputs'); root.mkdirs()
    def fa=root.resolve('input.fa'); fa.text='>c1\nACGTACGT\n>c2\nNNNN\n'
    def map=root.resolve('mapping'); map.mkdirs()
    map.resolve('coverage.txt').text='#rname\tstartpos\tendpos\tmeandepth\nc1\t1\t8\t20\nc2\t1\t4\t0\n'
    def assemblies=Channel.of(tuple([id:'A',sample:'S',taxid:7955],fa),tuple([id:"${'B'}",sample:'S',taxid:7955],fa))
    CLEAN(assemblies,Channel.value('/unused/test-db'),'contig')
    AUDIT(CLEAN.out.decontaminated,Channel.value('/unused/test-db'))
    PLOTS(AUDIT.out.assembly,AUDIT.out.taxonomy_report,Channel.empty(),'final',
        Channel.of(tuple([id:'A'],map),tuple([id:'B'],map)),true)
    REPORT(AUDIT.out.assembly,AUDIT.out.action_report,AUDIT.out.taxonomy_report,
        Channel.empty(),Channel.empty(),Channel.empty(),'final',true,false,PLOTS.out.metrics,AUDIT.out.adaptor_report)
    OVERVIEW(REPORT.out.tables)
    OVERVIEW.out.report.view { stage,p -> assert stage=='overview'; assert p.text.contains('FCS evidence'); 'PASS combined overview consumes current-run stage tables' }
    AUDIT.out.assembly.toList().view { rows ->
        assert rows.size()==2
        assert rows.every { it[1].text==fa.text }
        'PASS final audit preserves exact sequence bytes across two assemblies'
    }
    REPORT.out.report.view { stage,report ->
        assert report.text.contains('Final assembly audit')
        assert report.text.contains('No actionable FCS findings')
        assert report.text.contains('coloured by FCS-GX')
        assert !report.text.contains('## Optional BlobTools evidence')
        'PASS FCS-only report and plots without protein or nucleotide searches'
    }
}
