// Gather current-run evidence, never stale published directories.
process CONTAMINATION_COHORT_REPORT {
    conda 'conda-forge::python=3.11 conda-forge::r-base=4.3.3 conda-forge::r-ggplot2=3.5.2 conda-forge::r-scales=1.3.0'
    cpus 1
    memory '4 GB'
    time '1h'
    publishDir "${params.outdir}/${stage}/decontam", mode: 'copy'
    input:
    val stage
    val spec
    path artifacts, stageAs: 'inputs/input??/*'
    path helper
    path plotter
    output:
    tuple val(stage), path('contamination_report.md'), emit: report
    path 'report', emit: tables
    script:
    def plot = spec.fcs_metrics ? 'Rscript plot_contamination_cohort.R report/cohort_fcs_sequence_evidence.tsv report fcs' :
        (!spec.fcs_enabled && spec.blobtools_enabled && spec.blobtools ? 'Rscript plot_contamination_cohort.R report/cohort_sequence_evidence.tsv report' : '')
    """
    cat > report_spec.json <<'SPEC'
${groovy.json.JsonOutput.toJson(spec)}
SPEC
    python3 contamination_cohort_report.py --spec report_spec.json --out .
    ${plot}
    if [ -f report/cohort_fcs_segments.tsv ]; then
        Rscript plot_fcs_breakdown.R report/cohort_fcs_sequence_evidence.tsv report/cohort_fcs_segments.tsv report
    fi
    """
}

workflow CONTAMINATION_REPORTING {
    take:
    assemblies
    fcs_actions
    fcs_taxonomy
    fcs_contaminants
    blob_reviews
    blob_diagnostics
    stage
    fcs_enabled
    blobtools_enabled
    fcs_metrics
    adaptor_reports
    main:
    def fcs = fcs_actions.map { m,p -> tuple(m.id.toString(),p) }
        .join(fcs_taxonomy.map { m,p -> tuple(m.id.toString(),p) }, failOnMismatch:true)
        .map { id,a,t -> [assembly:id,actions:a.toString(),taxonomy:t.toString()] }.toList()
        .map { rows -> [rows:rows] }
        .combine(fcs_contaminants.toList().map { rows -> [rows:rows] })
        .map { a,c -> a.rows.collect { row ->
            def removed=c.rows.find { it[0].id.toString()==row.assembly }
            row + (removed ? [contaminants:removed[1].toString()] : [:]) } }
    def blob = blob_reviews.map { id,p -> tuple(id.toString(),p) }
        .join(blob_diagnostics.map { id,p -> tuple(id.toString(),p) }, failOnMismatch:true)
        .map { id,r,d -> [assembly:id.toString(),review:r.toString(),diagnostics:d.toString()] }.toList()
    def specification = assemblies.map { m,p -> m.id.toString() }.toList().map { ids -> [ids:ids] }
        .combine(fcs.map { records -> [records:records] }).combine(blob.map { records -> [records:records] })
        .map { ids,f,b -> [stage:stage,assemblies:ids.ids,fcs_enabled:fcs_enabled,blobtools_enabled:blobtools_enabled,fcs:f.records,blobtools:b.records] }
        .combine(fcs_metrics.toList().map { rows -> [metrics:rows.collect { id,p -> [assembly:id.toString(),metrics:p.toString()] }] })
        .combine(adaptor_reports.toList().map { rows -> [adaptor:rows.collect { m,p -> [assembly:m.id.toString(),report:p.toString()] }] })
        .map { s,m,a -> s + [fcs_metrics:m.metrics,adaptor_reports:a.adaptor] }
    CONTAMINATION_COHORT_REPORT(stage, specification,
        specification.map { s -> s.fcs.collectMany { [it.actions,it.taxonomy,it.contaminants].findAll { it }.collect { file(it) } } +
            s.blobtools.collectMany { [file(it.review),file(it.diagnostics)] } +
            s.fcs_metrics.collect { file(it.metrics) } + s.adaptor_reports.collect { file(it.report) } },
        ['contamination_cohort_report.py','fcs_category_evidence.py'].collect {
            file("${moduleDir}/../py_scripts/${it}",checkIfExists:true) },
        ['plot_contamination_cohort.R','fcs_plot_common.R','plot_fcs_breakdown.R'].collect {
            file("${moduleDir}/../r_scripts/${it}",checkIfExists:true) })
    emit:
    report = CONTAMINATION_COHORT_REPORT.out.report
    tables = CONTAMINATION_COHORT_REPORT.out.tables.map { p -> tuple(stage,p) }
}
