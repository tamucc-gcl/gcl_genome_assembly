/*
========================================================================================
    PANGENOME PRIVATE PLOTS MODULE
========================================================================================
    Repo location: modules/pangenome_private_plots.nf

    The private-sequence figure set.

    ONE TASK PER TAXID, TAKING BOTH FLAVOURS -- NOT FLAVOUR-PARALLEL
    ---------------------------------------------------------------
    Every other process in the private chain is flavour-parallel, but these figures are not,
    because the most important thing they have to convey is that the two arms DISAGREE:
    clipping removes 704,499,041 bp of which 698,360,436 -- 99.1% -- is private sequence. So
    the clip arm understates private content by 46%, and it also INVERTS the reference's
    apparent rank: highest of ten haplotypes on clip (15.08% of private bp), lowest of ten on
    full (8.08%), because the reference is the graph backbone and is never clipped.

    A figure drawn from one arm alone cannot show that, and two separate figures drawn one per
    arm leave the reader to do the comparison. The tables are small (a spectrum is ~10 rows
    per haplotype), so there is no cost to taking both.

    WHY THESE FIGURES ARE NOT IN PANGENOME_PLOTS
    -------------------------------------------
    PANGENOME_PLOTS is joined per-taxid with PANGENOME_GROWTH, which runs panacus on the CLIP
    GFA, so that process is structurally pinned to one flavour. Its two existing private
    figures were therefore drawn from the understating arm. They move here.

    Input : tuple(taxid, flavors, spectra, hap_privates, by_contigs, evidence_csvs, xtabs),
            reference id, script
    Output: figures / audit / versions

    Every input is optional at the script level -- a missing table skips its figure with a
    message rather than failing the task -- so a species where the private analysis was
    disabled still produces the rest of the report.
========================================================================================
*/

process PANGENOME_PRIVATE_PLOTS {
    tag "${taxid}"
    label 'pangenome_private_plots'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavors), path(spectra), path(hap_privates),
          path(by_contigs), path(evidence_csvs), path(xtabs)
    // CLIP ONLY, both of these. The SV catalog exists only on the clip graph, so a
    // full-arm tier spectrum would not be comparable to the SV panel above it -- and the
    // cactus authors recommend the clip graph for most applications.
    tuple val(ttaxid), path(tier_spectrum), path(sv_spectrum)
    tuple val(rtaxid), val(ref_id)
    path(script)

    output:
    // `${taxid}.*.png`, not `${taxid}.private_*.png`. The narrower glob collected only the
    // figures whose names happen to start with "private_", so size_by_sharing, copy_ratio and
    // evidence_by_chromosome were WRITTEN by the R and never published -- the R reported six
    // figures and Nextflow collected three, with no error either side.
    tuple val(taxid), path("${taxid}.*.png"), emit: figures, optional: true
    tuple val(taxid), path("${taxid}.private_figures_audit.tsv"), emit: audit
    path("versions.tsv"), emit: versions

    script:
    def minbp = params.pangenome_private_min_bp ?: 1000
    // Files arrive as flat lists in channel order; `flavors` is the parallel list of tags.
    // Pairing is done HERE from that list rather than by parsing filenames, because a
    // filename-derived flavour was already the source of one silent mismatch (the
    // `.private.fa` suffix that never stripped, which left every key wrong and the combine
    // matching nothing).
    def fl = (flavors instanceof List) ? flavors : [flavors]
    def asList = { it == null ? [] : ((it instanceof List) ? it : [it]) }
    def pick = { files, want ->
        def l = asList(files)
        def i = fl.findIndexOf { it == want }
        (i >= 0 && i < l.size()) ? l[i].name : 'NONE'
    }
    def kvs = []
    ['clip', 'full'].each { want ->
        kvs << "spectrum_${want}=${pick(spectra, want)}"
        kvs << "hap_${want}=${pick(hap_privates, want)}"
        kvs << "contig_${want}=${pick(by_contigs, want)}"
        kvs << "evidence_${want}=${pick(evidence_csvs, want)}"
        kvs << "xtab_${want}=${pick(xtabs, want)}"
    }
    def opts = kvs.join(' ')
    """
    set -euo pipefail

    echo "[PRIVATE_PLOTS ${taxid}] flavours present: ${fl.join(',')}" >&2

    Rscript ${script} \\
        ${taxid} . \\
        ${opts} \\
        tier_clip=${tier_spectrum} \\
        sv_clip=${sv_spectrum} \\
        ref='${ref_id}' \\
        min_bp=${minbp}

    # The audit records which figures were SKIPPED and why. A skipped figure is not a failure
    # -- a species with the private analysis disabled should still get the rest of the report
    # -- but a silently missing figure is how a broken input goes unnoticed, so the reasons
    # are surfaced in the log too.
    A=${taxid}.private_figures_audit.tsv
    if grep -q '^skipped_' "\$A"; then
        echo "[PRIVATE_PLOTS ${taxid}] figures skipped:" >&2
        grep '^skipped_' "\$A" | sed 's/^/  /' >&2
    fi
    # Same glob as the output declaration, so the count in the log matches what Nextflow
    # collects. They disagreed before -- the R reported six figures and this line said three.
    echo "[PRIVATE_PLOTS ${taxid}] \$(ls -1 ${taxid}.*.png 2>/dev/null | wc -l) figures" >&2

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tR\\t%s\\n' "${task.process}" "\$(Rscript --version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    printf 'metric\\tvalue\\nlabel\\t${taxid}\\n' > ${taxid}.private_figures_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
