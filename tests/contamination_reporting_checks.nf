nextflow.enable.dsl=2
include { CONTAMINATION_REPORTING as BOTH } from '__REPO__/workflows/contamination_reporting.nf'
include { CONTAMINATION_REPORTING as FCS_ONLY } from '__REPO__/workflows/contamination_reporting.nf'
include { CONTAMINATION_REPORTING as DISABLED } from '__REPO__/workflows/contamination_reporting.nf'
params.outdir='/tmp/contamination-report-check/results'
workflow {
    def root=file('/tmp/contamination-report-check/inputs')
    root.mkdirs()
    def action=root.resolve('A.fcs_gx_report.txt')
    action.text='## metadata\n#seq_id\tstart_pos\tend_pos\taction\nc1\t1\t10\tEXCLUDE\n'
    def taxonomy=root.resolve('A.taxonomy.rpt'); taxonomy.text='taxonomy\n'
    def contaminants=root.resolve('A.contaminants.fasta'); contaminants.text='>c1\nACGT\n'
    def fasta=root.resolve('A.fasta'); fasta.text='>c2\nACGT\n'
    def review=root.resolve('review'); review.mkdirs()
    review.resolve('evidence.json').text=groovy.json.JsonOutput.toJson([assembly:'A',whole_contig_candidates:0,
        protein_searched_bp:4,megablast_searched_bp:4,nucleotide_status:'not-run'])
    review.resolve('review.tsv').text='action_id\tsupport\nA1\tMEGABLAST\n'
    review.resolve('sequence_evidence.tsv').text='sequence\tlength_bp\tgc_fraction\tmean_depth\tplot_taxon\nc2\t4\t0.5\t20\tChordata\n'
    def diagnostics=root.resolve('diagnostics'); diagnostics.mkdirs(); diagnostics.resolve('blob.png').text='fixture'
    def fcsMetrics=root.resolve('A.fcs_metrics.tsv')
    fcsMetrics.text='sequence\tlength_bp\tgc_fraction\tmean_depth\tplot_taxon\tclassification_source\nc2\t4\t0.5\t20\tanml:fishes\tFCS-GX\n'
    def adaptor=root.resolve('A.fcs_adaptor_report.txt')
    adaptor.text='#accession\tlength\taction\trange\tname\n'
    BOTH(Channel.value(tuple([id:'A'],fasta)),Channel.value(tuple([id:'A'],action)),
        Channel.value(tuple([id:'A'],taxonomy)),Channel.value(tuple([id:'A'],contaminants)),
        Channel.value(tuple('A',review)),Channel.value(tuple('A',diagnostics)), 'contig',true,true,
        Channel.value(tuple('A',fcsMetrics)),Channel.value(tuple([id:'A'],adaptor)))
    FCS_ONLY(Channel.value(tuple([id:'A'],fasta)),Channel.value(tuple([id:'A'],action)),
        Channel.value(tuple([id:'A'],taxonomy)),Channel.value(tuple([id:'A'],contaminants)),
        Channel.empty(),Channel.empty(),'scaffold',true,false,
        Channel.value(tuple('A',fcsMetrics)),Channel.value(tuple([id:'A'],adaptor)))
    DISABLED(Channel.value(tuple([id:'A'],fasta)),Channel.empty(),Channel.empty(),Channel.empty(),
        Channel.empty(),Channel.empty(),'disabled',false,false,Channel.empty(),Channel.empty())
    BOTH.out.report.view { stage,report ->
        assert report.text.contains('[GX actions](fcs/A.fcs_gx_report.txt)')
        assert report.text.contains('[Decision sheet](blobtools/A/review/review.tsv)')
        'PASS combined report: both evidence sources'
    }
    FCS_ONLY.out.report.view { stage,report ->
        assert report.text.contains('Optional BlobTools searches: **not run**')
        assert report.text.contains('coloured by FCS-GX')
        'PASS combined report: FCS only'
    }
    DISABLED.out.report.view { stage,report ->
        assert report.text.contains('FCS: **not run**')
        'PASS combined report: disabled sources'
    }
}
