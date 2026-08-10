/*
========================================================================================
    PANGENOME REPORT MODULE  (workstream F)
========================================================================================
    Repo location: modules/pangenome_report.nf

    Renders a self-contained pangenome report section (markdown) + a machine-readable
    stats JSON from the pangenome stats tables, via r_scripts/pangenome_report.R. The
    markdown is written in the same style as generate_summary_report.R so the main report
    can append it as a child section (see the integration note below).

    Always on when a pangenome is built (gated by params.pangenome_report). Reuses the
    report R stack via the 'summarize_assembly' label. Inputs that are absent are passed
    as NO_* sentinels and skipped by the R script.

    Main-report integration (apply in main.nf / reporting.nf / generate_summary_report.R):
      - generate_summary_report.R: add `--pangenome_report` arg; before writing, if it is
        not a sentinel/empty, `md <- c(md, "", readLines(pangenome_report))`.
      - SUMMARY_REPORT: add `path(pangenome_report)` input + `--pangenome_report ...` arg.
      - REPORTING: add a `ch_pangenome_report` take and pass it through.
      - main.nf: pass PANGENOME.out.report (or a NO_PANGENOME sentinel) into REPORTING.

    Input : tuple(taxid, qc_metrics, growth_fit, variant_summary, graph_stats), report_script
    Output: report (pangenome_report.md) / json (pangenome_stats.json) / versions
========================================================================================
*/

process PANGENOME_REPORT {
    tag "taxid_${taxid}"
    label 'summarize_assembly'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(qc_metrics), path(growth_fit), path(variant_summary), path(graph_stats), path(popstruct), path(progressive)
    path(report_script)

    output:
    tuple val(taxid), path("pangenome_report.md"),   emit: report
    tuple val(taxid), path("pangenome_stats.json"),  emit: json
    path("versions.tsv"),                            emit: versions

    script:
    """
    Rscript ${report_script} \\
        --qc_metrics ${qc_metrics} \\
        --growth_fit ${growth_fit} \\
        --variant_summary ${variant_summary} \\
        --graph_stats ${graph_stats} \\
        --popstruct ${popstruct} \\
        --progressive ${progressive} \\
        --species ${taxid} \\
        --output pangenome_report.md \\
        --json pangenome_stats.json

    printf 'process\\ttool\\tversion\\n%s\\tRscript\\t%s\\n' "${task.process}" "\$(Rscript --version 2>&1 | awk '{print \$NF}')" > versions.tsv
    """

    stub:
    """
    printf '## Pangenome — %s\\n' "${taxid}" > pangenome_report.md
    printf '{}\\n' > pangenome_stats.json
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
