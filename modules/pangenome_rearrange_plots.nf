/*
========================================================================================
    PANGENOME REARRANGE PLOTS MODULE
========================================================================================
    Repo location: modules/pangenome_rearrange_plots.nf

    Layer 1 figures: inversions, duplications, rearrangement candidates.

    FLAVOUR-PARALLEL, unlike PANGENOME_PRIVATE_PLOTS
    -----------------------------------------------
    The private figures take both arms in one task because their whole point is that the arms
    disagree about how much private sequence exists. Rearrangement is different: the full arm
    is simply the correct one to read, because clipping cuts paths into subpaths (556 vs 394 on
    chr10) and a rearrangement straddling a boundary is lost to path projection entirely. The
    clip arm is run for comparison, not because its numbers are wanted, so the two are drawn
    separately and labelled by flavour rather than overlaid.

    WHY THESE FIGURES CARRY THE REARRANGEMENT ARGUMENT
    -------------------------------------------------
    Inversions in this graph are overwhelmingly NOT bubble-representable -- 29 path-explicit
    alleles and 411 alignment-rescued, against ~21.5 Mb of inverted sequence on chr10 alone
    from path projection -- and `self.cov > 1` is the ONLY duplication signal anywhere, since
    AT traversals found 27 node re-visits in 3,268,312 SV alleles. For rearrangement these are
    not a second opinion on the variant catalog; they are the entire evidence base.

    THE CANDIDATE FIGURE DEPENDS ON THE CARRIER SPLIT
    ------------------------------------------------
    It plots span against n_chrom_INDIVIDUALS, which only exists in candidate tables written
    by the current rearrange_from_untangle.py. The R checks for that column and skips the
    figure with an explicit "rerun REARRANGE" message rather than silently plotting the wrong
    axis, because the pre-split table's n_carriers pooled 74 unplaced projections with 6 real
    chromosome-scale carriers at the chr10 locus.

    Input : tuple(taxid, flavor, candidates, inversions, duplications, orientation, audit),
            script
    Output: figures / audit / versions

    Inputs are optional at the script level: a missing table skips its figure with a message
    rather than failing the task, and every skip is recorded in the audit and echoed into the
    log so a quietly absent figure is still visible.
========================================================================================
*/

process PANGENOME_REARRANGE_PLOTS {
    tag "${taxid}:${flavor}"
    label 'pangenome_rearrange_plots'

    publishDir "${params.outdir}/pangenome/${taxid}/rearrange", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), path(candidates), path(inversions),
          path(duplications), path(orientation), path(untangle_audit)
    path(script)

    output:
    tuple val(taxid), val(flavor), path("${taxid}.${flavor}.*.png"),
        emit: figures, optional: true
    tuple val(taxid), val(flavor), path("${taxid}.${flavor}.rearrange_figures_audit.tsv"),
        emit: audit
    path("versions.tsv"), emit: versions

    script:
    // min_span is the same floor the untangle run used, so the figures describe the same
    // population as the tables rather than a differently-filtered one.
    def minspan = params.pangenome_inv_min_bp ?: 1000
    """
    set -euo pipefail

    Rscript ${script} \\
        ${taxid}.${flavor} . \\
        candidates=${candidates} \\
        inversions=${inversions} \\
        duplications=${duplications} \\
        orientation=${orientation} \\
        audit=${untangle_audit} \\
        min_span=${minspan}

    A=${taxid}.${flavor}.rearrange_figures_audit.tsv
    if grep -q '^skipped_' "\$A"; then
        echo "[REARRANGE_PLOTS ${taxid}:${flavor}] figures skipped:" >&2
        grep '^skipped_' "\$A" | sed 's/^/  /' >&2
    fi

    # The candidate figure needs n_chrom_individuals, which only the current
    # rearrange_from_untangle.py emits. Surface that specifically: it means REARRANGE ran with
    # a stale script, and the pre-split table pooled 74 unplaced projections with 6 real
    # chromosome-scale carriers at the chr10 locus.
    if grep -q 'predates the carrier split' "\$A"; then
        echo "[REARRANGE_PLOTS ${taxid}:${flavor}] WARNING: candidate table predates the" >&2
        echo "  carrier split. Update py_scripts/rearrange_from_untangle.py and rerun" >&2
        echo "  PANGENOME_REARRANGE." >&2
    fi

    echo "[REARRANGE_PLOTS ${taxid}:${flavor}] \$(ls -1 ${taxid}.${flavor}.*.png 2>/dev/null | wc -l) figures" >&2

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tR\\t%s\\n' "${task.process}" "\$(Rscript --version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    printf 'metric\\tvalue\\nlabel\\t${taxid}.${flavor}\\n' \\
      > ${taxid}.${flavor}.rearrange_figures_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
