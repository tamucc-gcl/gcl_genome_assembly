/*
========================================================================================
    PANGENOME REPORT FIGURES MODULE  (workstream D)
========================================================================================
    Repo location: modules/pangenome_plots.nf

    Renders the pangenome report figures from the panacus coverage histogram and the
    variant catalog, via r_scripts/pangenome_plots.R. Growth/core curves, the Heaps'-law
    fit, and the confidence band are computed from the coverage histogram (rarefaction
    formulas); the SV size spectrum and variant-class bar come from the catalog.

    Reuses the pipeline R stack (ggplot2) via the 'pairwise_alignment' label.

    Input : tuple(taxid, hist_tsv, sv_sizes_tsv, variant_summary_tsv), plots_script
    Output: figures (PDFs) / growth_fit (machine-readable gamma / open-closed / sizes)
========================================================================================
*/

process PANGENOME_PLOTS {
    tag "${taxid}"
    label 'pairwise_alignment'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    // hap_private is GONE: this process is pinned to the clip arm by its join with
    // PANGENOME_GROWTH (panacus on the clip GFA), and the private figures need both arms.
    // They live in PANGENOME_PRIVATE_PLOTS now.
    //
    // footprint / length_class / ref_fai may be NO_FILE -- the R skips a figure whose table is
    // missing rather than failing, so a species with CLASSIFY disabled still gets the growth
    // and coverage figures.
    tuple val(taxid), path(hist), path(sv_sizes), path(variant_summary),
          path(footprint), path(length_class), path(ref_fai)
    path(plots_script)

    output:
    tuple val(taxid), path("${taxid}.*.png"),          emit: figures, optional: true
    tuple val(taxid), path("${taxid}.growth_fit.tsv"), emit: growth_fit

    script:
    """
    Rscript ${plots_script} ${hist} ${sv_sizes} ${variant_summary} ${taxid} . \\
        footprint=${footprint} \\
        length_class=${length_class} \\
        fai=${ref_fai} \\
        core=${params.pangenome_tier_core} \\
        softcore=${params.pangenome_tier_softcore} \\
        shell=${params.pangenome_tier_shell}

    # A footprint spanning ONE chromosome means the coordinate frames are still being mixed --
    # the bug that gave SUBST 89.5 Mb where the correct value is 420.2 Mb. The R reports the
    # count; surface it here so it appears in the task log rather than only in stderr.
    if [ -s ${footprint} ]; then
        n=\$(awk -F'\\t' '!/^#/ && \$1!="chrom"{print \$1}' ${footprint} | sort -u | wc -l)
        echo "[PLOTS ${taxid}] reference footprint spans \$n chromosome(s)" >&2
        if [ "\${n:-0}" -le 1 ]; then
            echo "[PLOTS ${taxid}] WARNING: a single chromosome means the footprint frames" >&2
            echo "  are mixed. Check classify_variants.py keys ref_iv by (class, chrom)." >&2
        fi
    fi
    """

    stub:
    """
    : > ${taxid}.growth_curves.png
    printf 'metric\\tvalue\\n' > ${taxid}.growth_fit.tsv
    """
}
