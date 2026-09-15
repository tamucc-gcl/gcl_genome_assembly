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

    Input : tuple(taxid, qc_metrics, growth_fit, variant_summary, graph_stats, popstruct,
                  progressive, manifest, hap_private), report_script
    Output: report (pangenome_report.md) / json (pangenome_stats.json) / versions
========================================================================================
*/

process PANGENOME_REPORT {
    tag "taxid_${taxid}"
    label 'summarize_assembly'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    // The last six feed the "Which view says what" matrix. All tolerate a NO_* placeholder:
    // the R degrades a cell to "—" rather than failing, so a species built without CLASSIFY
    // or without the private analysis still gets a report, with the matrix showing which
    // views were not built.
    tuple val(taxid), path(qc_metrics), path(growth_fit), path(variant_summary),
          path(graph_stats), path(popstruct), path(progressive), path(manifest),
          path(hap_private),
          path(audit_parent), path(audit_fine), path(hap_private_full),
          path(priv_figures_audit), path(rearr_audit_clip), path(rearr_audit_full),
          path(input_cov_audit)
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
        --manifest ${manifest} \\
        --hap_private ${hap_private} \\
        --audit_parent ${audit_parent} \\
        --audit_fine ${audit_fine} \\
        --hap_private_full ${hap_private_full} \\
        --priv_figures_audit ${priv_figures_audit} \\
        --rearr_audit_clip ${rearr_audit_clip} \\
        --rearr_audit_full ${rearr_audit_full} \\
        --input_cov_audit ${input_cov_audit} \\
        --species ${taxid} \\
        --output pangenome_report.md \\
        --json pangenome_stats.json

    # The matrix is only worth reading if the topology row reflects what the classifier
    # actually did. Surface it in the task log too, because a report claiming topological
    # classes on a decomposed VCF would be wrong in a way nobody would catch by reading it.
    if [ -s ${audit_fine} ]; then
        te=\$(awk -F'\\t' '\$1=="topology_enabled"{print \$2}' ${audit_fine})
        echo "[REPORT ${taxid}] fine-tier topology_enabled=\${te:-unknown}" >&2
        if [ "\${te:-}" = "True" ]; then
            echo "[REPORT ${taxid}] WARNING: the fine tier reports topology ENABLED. It is" >&2
            echo "  vcfwave-decomposed, so AT is inherited from the parent record while" >&2
            echo "  REF/ALT are rewritten -- allele i no longer matches traversal i+1." >&2
            echo "  Check the decomposition markers in classify_variants.py." >&2
        fi
    fi

    printf 'process\\ttool\\tversion\\n%s\\tRscript\\t%s\\n' "${task.process}" "\$(Rscript --version 2>&1 | awk '{print \$NF}')" > versions.tsv
    """

    stub:
    """
    printf '## Pangenome — %s\\n' "${taxid}" > pangenome_report.md
    printf '{}\\n' > pangenome_stats.json
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
