/*
========================================================================================
    PANGENOME 2D LAYOUT MODULE  (workstream D)
========================================================================================
    Repo location: modules/pangenome_2d_viz.nf

    2D graph layout (odgi sort -O -> odgi layout -> odgi draw) for ONE chromosome graph,
    as a complement to the 1D odgi viz that cactus already produces. The subworkflow
    scatters this over the clip per-chromosome .og graphs, so the (expensive, ~25 min each)
    path-guided SGD layouts run as parallel tasks rather than one long serial loop.

    odgi layout/draw build an index that requires OPTIMIZED (contiguous) node IDs; cactus's
    .og is not, so each chromosome is optimized on a throwaway copy first (odgi sort -O).

    Best-effort: errorStrategy 'ignore' so a slow/failed chromosome (e.g. a per-task time
    limit on a very large graph) never kills the pipeline; the PNGs that do complete are
    still published. Gated by params.pangenome_2d_viz (default on). Layout cost grows with
    graph size and haplotype count, so disable it for very large production runs if it
    dominates runtime.

    Input : tuple(taxid, og)   (one clip per-chromosome .og)
    Output: png (this chromosome's 2D PNG) / versions
========================================================================================
*/

process PANGENOME_2D_VIZ {
    tag "${taxid}:${og.simpleName}"
    label 'cactus_tools'

    // per-chromosome, best-effort: one slow/oversized chromosome must not fail the run
    errorStrategy 'ignore'
    time '4h'

    publishDir "${params.outdir}/pangenome/${taxid}/layout2d", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(og)

    output:
    tuple val(taxid), path("${og.simpleName}.2D.png"), emit: png,      optional: true
    path("versions.tsv"),                              emit: versions, optional: true

    script:
    def base = og.simpleName
    """
    set -euo pipefail
    export HOME="\$PWD"

    # optimize node IDs (required by layout/draw index), then 2D layout + render.
    # -P prints layout progress (verbose in .command.err, but lets you monitor the ~25-min SGD)
    odgi sort   -i ${og} -o ${base}.sorted.og -O
    odgi layout -i ${base}.sorted.og -o ${base}.lay -t ${task.cpus} -P
    odgi draw   -i ${base}.sorted.og -c ${base}.lay -p ${base}.2D.png -H 1500

    rm -f ${base}.sorted.og ${base}.lay

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\todgi\\t%s\\n' "${task.process}" "\$(odgi version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    : > ${og.simpleName}.2D.png
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
