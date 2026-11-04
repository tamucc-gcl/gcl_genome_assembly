include { REPORT_PAINTING } from '../modules/report_painting.nf'
include { MAIN_REPORT } from '../modules/main_report.nf'
include { ANNOTATION_AUDIT } from '../modules/annotation_audit.nf'
include { graphReportPath; graphReviewRows } from '../functions/report_artifacts.nf'

workflow REPORT_PAINTINGS {
    take:
    eligibility
    main:
    blocks=Channel.empty();versions=Channel.empty()
    if(params.run_report && params.report_chromosome_painting) {
        jobs=eligibility.flatMap { manifest ->
            def data=new groovy.json.JsonSlurper().parse(manifest)
            def refs=data.linear_references.collectEntries { [(it.taxid.toString()):it] }
            data.assemblies.findAll { it.completion=='finalized' && refs[it.taxid.toString()] }.collect { row ->
                def ref=refs[row.taxid.toString()]
                tuple([id:row.id,sample:row.sample,taxid:row.taxid.toString(),haplotype:row.haplotype],ref.id,
                    file(ref.fasta,checkIfExists:true),file(ref.fai,checkIfExists:true),
                    file(row.fasta,checkIfExists:true),file(row.fai,checkIfExists:true))
            }
        }
        REPORT_PAINTING(jobs,file("${projectDir}/py_scripts/report_painting.py",checkIfExists:true))
        blocks=REPORT_PAINTING.out.blocks
        versions=REPORT_PAINTING.out.versions
    }
    emit:
    blocks=blocks
    versions=versions
}

workflow MAIN_REPORTING {
    take:
    eligibility
    expected
    name_maps
    paintings
    qc_metrics
    decisions
    verifications
    graph_reports
    detailed_report
    dotplots
    ribbons
    contacts
    organelle_maps
    organelle_stats
    annotations
    recovery_status
    qc_plots
    finishing_stats
    qc_html
    final_snails
    input_capabilities
    chimera_registry
    versions
    main:
    def artifact={ type,id,path -> [type:type,id:id.toString(),path:path.toString()] }
    def annotation_helper=file("${projectDir}/py_scripts/annotation_counts.py",checkIfExists:true)
    ANNOTATION_AUDIT(annotations,annotation_helper)
    registry=graph_reports.map { taxid,md ->
        def review=md.name in ['source_sequences.tsv','source_audit.json']
        [type:md.name=='sharing_summary.tsv' ? 'pangenome_summary' : md.name=='source_sequences.tsv' ? 'graph_review_rows' : md.name=='source_audit.json' ? 'graph_review_audit' : 'pangenome',id:taxid.toString(),
         path:review ? "pangenome/${taxid}/representation/${md.name}".toString() : md.name=='sharing_summary.tsv' ? "pangenome/${taxid}/sharing/sharing_summary.tsv".toString() : graphReportPath(taxid,md.name),
         content:md.name=='source_sequences.tsv' ? graphReviewRows(md) : review || md.name=='sharing_summary.tsv' ? md.text : '',
         note:md.name=='source_report.md' ? (md.readLines().find { it.startsWith('Status:') } ?: '') : '']
    }.mix(detailed_report.map { md -> artifact('assembly_details','run',md.name) })
      .mix(dotplots.map { a,b,p -> artifact('dotplot',"${a} / ${b}","pairwise_alignments/${p.name}") })
      .mix(ribbons.map { a,b,p -> artifact('ribbon',"${a} / ${b}","pairwise_alignments/${p.name}") })
      .mix(contacts.flatMap { m,stage,files -> (files instanceof List ? files : [files]).collect { p ->
          artifact('contact_map',m.id,"contact_maps/${p.name}") } })
      .mix(organelle_maps.map { m,p -> artifact('organelle',m.sample,"organelle/${p.name}") })
      .mix(organelle_stats.map { m,p -> [type:'organelle_stats',id:m.sample.toString(),path:"organelle/${p.name}".toString(),content:p.text] })
      .mix(annotations.map { m,p -> artifact('organelle_annotation',m.sample,"organelle/${p.name}") })
      .mix(ANNOTATION_AUDIT.out.counts.map { m,p -> [type:'annotation_counts',id:m.sample.toString(),path:"reports/organelle_annotations/${p.name}".toString(),content:p.text] })
      .mix(ANNOTATION_AUDIT.out.loci.map { m,p -> artifact('annotation_loci',m.sample,"reports/organelle_annotations/${p.name}") })
      .mix(recovery_status.map { m,org,p -> [type:'organelle_recovery',id:m.sample.toString(),organelle:org.toString(),path:"organelle/${p.name}".toString(),content:p.text] })
      .mix(qc_plots.flatMap { files -> (files instanceof List ? files : [files]).collect { p -> artifact('stage_qc',p.name,"qc/assembly/${p.name}") } })
      .mix(finishing_stats.filter { p -> !p.name.startsWith('NO_') }.map { p -> artifact('telomere_extension','run',"assembly/scaffold/teloclip/${p.name}") })
      .mix(qc_html.map { p -> artifact('interactive_qc','run',"reports/${p.name}") })
      .mix(final_snails.map { id,stage,p -> artifact('snail',id,"snail_plots/${p.name}") })
      .mix(chimera_registry.filter { p -> !p.name.startsWith('NO_') }.map { p ->
          [type:'chimera_registry',id:'run',path:'assembly/chimeras/assembly_registry.tsv',content:p.text] })
      .mix(versions.flatMap { files -> (files instanceof Collection ? files : [files]).findAll { !it.name.startsWith('NO_') }.collect { p ->
          artifact('provenance',p.name,"provenance/${p.name}") } })
      .toList()
    def fallback=file("${projectDir}/assets/NO_REPORT_INPUT",checkIfExists:true)
    context=input_capabilities.map { samples -> [pipeline:workflow.manifest.name ?: 'gcl_genome_assembly',
        version:workflow.manifest.version,revision:workflow.revision,commit:workflow.commitId,
        run_name:workflow.runName.toString(),profile:workflow.profile.toString(),
        nextflow:nextflow.version.toString(),start:workflow.start.toString(),
        qc_mode:(params.qc_mode ?: 'all_stages').toString(),chimera_enabled:params.chimera_detect,
        painting_enabled:params.report_chromosome_painting,top_scaffolds:params.report_top_scaffolds,
        contract_enabled:true,samples:samples,contact_maps_enabled:params.run_final_contact_maps != false,
        contact_resolutions:(params.hic_plot_resolutions ?: '1000000,500000,250000,100000').toString(),
        comparisons_enabled:params.run_post_assembly != false && params.run_pairwise_alignments,
        comparison_mode:(params.pairwise_alignment_mode ?: 'within_sample').toString(),
        graph_enabled:params.run_pangenome,graph_validate_only:params.pangenome_validate_only] }
    MAIN_REPORT(eligibility,expected,name_maps.map { m,p -> p }.toList(),
        paintings.flatMap { m,ref,tsv,json -> [tsv,json] }.toList(),
        qc_metrics.ifEmpty(fallback),decisions,verifications,registry,context,versions,
        file("${projectDir}/r_scripts/generate_main_report.R",checkIfExists:true))
    emit:
    report=MAIN_REPORT.out.report
}
