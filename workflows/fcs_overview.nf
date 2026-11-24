process FCS_EVIDENCE_OVERVIEW {
    conda 'conda-forge::python=3.11 conda-forge::r-base=4.3.3 conda-forge::r-ggplot2=3.5.2 conda-forge::r-scales=1.3.0'
    cpus 1
    memory '4 GB'
    time '1h'
    publishDir "${params.outdir}/overview/decontam", mode:'copy',saveAs:{ name -> name.startsWith('overview/') ? name.substring(9) : name }
    input:
    path directories, stageAs:'inputs/input??/*'
    path helper
    path renderers
    output:
    tuple val('overview'),path('overview/contamination_report.md'),emit:report
    path 'overview/*',emit:artifacts
    script:
    def q={ p -> "'"+p.toString().replace("'","'\"'\"'")+"'" }
    def inputs=(directories instanceof Collection ? directories : [directories]).collect { q(it) }.join(' ')
    """
    python3 fcs_overview.py ${inputs}
    if [ -f overview/cohort_fcs_sequence_evidence.tsv ]; then
        Rscript plot_fcs_overview.R overview/cohort_fcs_sequence_evidence.tsv overview/cohort_fcs_composition
    fi
    """
}
workflow FCS_OVERVIEW {
    take:
    tables
    main:
    FCS_EVIDENCE_OVERVIEW(tables.map { stage,p -> p }.toList(),
        file("${moduleDir}/../py_scripts/fcs_overview.py",checkIfExists:true),
        ['fcs_plot_common.R','plot_fcs_overview.R'].collect { file("${moduleDir}/../r_scripts/${it}",checkIfExists:true) })
    emit:
    report=FCS_EVIDENCE_OVERVIEW.out.report
}
