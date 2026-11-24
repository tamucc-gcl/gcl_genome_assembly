nextflow.enable.dsl=2
include { DECONTAMINATE_ASSEMBLY as CLEAN } from '../../workflows/decontaminate_assembly.nf'
include { FCS_FINAL_SCREEN as AUDIT } from '../../workflows/fcs_final_screen.nf'
include { FCS_DIAGNOSTICS as CONTIG_PLOTS; FCS_DIAGNOSTICS as FINAL_PLOTS } from '../../workflows/fcs_diagnostics.nf'
include { CONTAMINATION_REPORTING as CONTIG_REPORT; CONTAMINATION_REPORTING as FINAL_REPORT } from '../../workflows/contamination_reporting.nf'

workflow {
    if(params.run_contamination_evidence) error 'This bounded test must not enable BlobTools searches'
    def meta=[id:params.assembly,sample:params.sample,taxid:params.taxid,hifi:true]
    def source=Channel.value(tuple(meta,file(params.source_fasta,checkIfExists:true)))
    def finalAssembly=Channel.value(tuple(meta,file(params.final_fasta,checkIfExists:true)))
    def reads=Channel.value(tuple(meta,file(params.hifi_fastq,checkIfExists:true)))
    def database=Channel.value(params.gxdb_dir)
    CLEAN(source,database,'contig')
    CONTIG_PLOTS(CLEAN.out.screened_input,CLEAN.out.taxonomy_report,reads,'contig',Channel.empty(),false)
    CONTIG_REPORT(CLEAN.out.decontaminated,CLEAN.out.action_report,CLEAN.out.taxonomy_report,
        CLEAN.out.contaminants,Channel.empty(),Channel.empty(),'contig',true,false,
        CONTIG_PLOTS.out.metrics,CLEAN.out.adaptor_report)
    // Audit the actual retained finalized assembly, not the newly cleaned contigs.
    AUDIT(finalAssembly,database)
    FINAL_PLOTS(AUDIT.out.assembly,AUDIT.out.taxonomy_report,reads,'final',Channel.empty(),false)
    FINAL_REPORT(AUDIT.out.assembly,AUDIT.out.action_report,AUDIT.out.taxonomy_report,
        Channel.empty(),Channel.empty(),Channel.empty(),'final',true,false,
        FINAL_PLOTS.out.metrics,AUDIT.out.adaptor_report)
}
