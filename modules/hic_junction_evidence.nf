/*
========================================================================================
    HI-C JUNCTION EVIDENCE
========================================================================================
    Repo location: modules/hic_junction_evidence.nf

    Per assembly: does Hi-C say two of its scaffolds are one chromosome the scaffolder
    failed to join? A WITHIN-assembly measurement, so it is independent of the
    harmonization frame, of the reference choice, and of every other assembly -- which is
    what lets it break ties the cross-assembly alignment consensus cannot.

    Runs on PRE-finalize pairs, in the same coordinates as the assembly HARMONIZE is about
    to harmonize. That is the whole reason HIC_PREFINAL_PAIRS exists upstream of
    HARMONIZE_SCAFFOLDS.

    Candidate scaffolds come from this assembly's own fai by size threshold, deliberately
    NOT from the drop-off caller -- depending on the frame would be circular, and
    over-inclusion only widens the null.

    Input : tuple(meta, pairs_gz, fai) ; path(script)
    Output: evidence (TSV) / versions
========================================================================================
*/

nextflow.enable.dsl = 2

process HIC_JUNCTION_EVIDENCE {
    tag "${meta.id}"
    label 'hic_junction_evidence'

    publishDir "${params.outdir}/assembly/harmonization/hic_evidence",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.tsv' ? null : fn }

    input:
    tuple val(meta), path(pairs_gz), path(fai)
    path(script)

    output:
    tuple val(meta), path("${meta.id}.junction_evidence.tsv"), emit: evidence
    path("versions.tsv"), emit: versions

    script:
    def min_scaf = params.hic_junction_min_scaffold_bp ?: (params.finalize_min_scaffold_bp ?: 1000000)
    def wfrac    = (params.hic_junction_window_frac   != null) ? params.hic_junction_window_frac   : 0.25
    def wmax     = (params.hic_junction_window_max_bp != null) ? params.hic_junction_window_max_bp : 10000000
    def minenr   = (params.hic_junction_min_enrichment != null) ? params.hic_junction_min_enrichment : 2.0
    def minz     = (params.hic_junction_min_z != null) ? params.hic_junction_min_z : 5.0
    def mintot   = (params.hic_junction_min_contacts != null) ? params.hic_junction_min_contacts : 100
    def maxcand  = (params.hic_junction_max_candidates != null) ? params.hic_junction_max_candidates : 1000
    """
    set -euo pipefail

    python3 ${script} \\
        --pairs ${pairs_gz} \\
        --fai ${fai} \\
        --assembly-id ${meta.id} \\
        --out ${meta.id}.junction_evidence.tsv \\
        --min-scaffold-bp ${min_scaf} \\
        --window-frac ${wfrac} \\
        --window-max-bp ${wmax} \\
        --min-corner-enrichment ${minenr} \\
        --min-z ${minz} \\
        --min-total-contacts ${mintot} \\
        --max-candidates ${maxcand}

    PY=\$(python3 --version 2>&1 | awk '{print \$2}')
    {
        printf 'process\\ttool\\tversion\\n'
        printf '%s\\tpython\\t%s\\n' "${task.process}" "\${PY}"
    } > versions.tsv
    """

    stub:
    """
    printf 'assembly\\tscaf_a\\tscaf_b\\tlen_a\\tlen_b\\twindow_a\\twindow_b\\ttotal_contacts\\tbest_corner\\tcorner_contacts\\texpected_contacts\\tcorner_enrichment\\tpoisson_z\\tcall\\torientation\\n' \\
        > ${meta.id}.junction_evidence.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
