/*
========================================================================================
    PANGENOME 2D LAYOUT MODULE  (workstream D)
========================================================================================
    Repo location: modules/pangenome_2d_viz.nf

    Per-chromosome 2D graph layout (odgi sort -> odgi layout -> odgi draw) as a complement
    to the 1D odgi viz that cactus already produces. Operates on the clip per-chromosome .og
    graphs. odgi layout/draw need a sorted graph, so each chromosome is sorted on a throwaway
    copy first (the sort is a viz concern, kept out of the build step). Best-effort: any
    per-chromosome failure is skipped so it never breaks the pipeline (odgi is in the cactus
    image; the 'png' output is optional).

    Gated by params.pangenome_2d_viz (default on). Layout is iterative (path-guided SGD),
    so cost grows with graph size — disable for very large runs if it dominates runtime.

    Input : tuple(taxid, [chrom_ogs])   (CACTUS_PANGENOME.out.chrom_og; clip + full .og)
    Output: png (per-chromosome 2D PNGs) / versions
========================================================================================
*/

process PANGENOME_2D_VIZ {
    tag "${taxid}"
    label 'cactus_tools'

    publishDir "${params.outdir}/pangenome/${taxid}/layout2d", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(chrom_ogs)

    output:
    tuple val(taxid), path("*.2D.png"), emit: png, optional: true
    path("versions.tsv"),               emit: versions

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"

    for og in ${chrom_ogs}; do
        case "\$og" in *.full.og) continue;; esac        # clip per-chr graphs only
        base=\$(basename "\$og" .og)
        # odgi layout/draw need a sorted (optimized) graph; cactus's .og is not sorted for
        # 2D layout, so sort a throwaway copy first (path-guided SGD + groom + topological).
        odgi sort   -i "\$og" -o "\$base.sorted.og" -p Ygs -t ${task.cpus} -P || { echo "sort failed: \$base"   >&2; continue; }
        odgi layout -i "\$base.sorted.og" -o "\$base.lay" -t ${task.cpus} -P  || { echo "layout failed: \$base" >&2; rm -f "\$base.sorted.og"; continue; }
        odgi draw   -i "\$base.sorted.og" -c "\$base.lay" -p "\$base.2D.png" -H 1500 || echo "draw failed: \$base" >&2
        rm -f "\$base.sorted.og" "\$base.lay"
    done

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\todgi\\t%s\\n' "${task.process}" "\$(odgi version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    : > chr1_1.2D.png
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
