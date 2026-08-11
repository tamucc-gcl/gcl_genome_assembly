/*
========================================================================================
    PANGENOME PROGRESSIVE GROWTH MODULE  (workstream H — opt-in)
========================================================================================
    Repo location: modules/pangenome_progressive.nf

    EMPIRICAL pangenome growth by incremental construction: start from the reference and
    add one assembly at a time, recording the graph size (nodes / edges / bp) after each
    addition. This complements the analytic panacus growth (workstream E), which is the
    exact expected growth from the FINISHED graph — here we measure the graph as it is
    actually built, in the deterministic reference-first order.

    Built on minigraph (the SV-graph construction backbone of Minigraph-Cactus), which
    supports incremental extension (`minigraph -cxggs existing.gfa new.fa`) and is cheap.
    A full cactus rebuild per step would be the base-level analogue but is prohibitively
    expensive; this is minigraph-level (SV) growth. Opt-in via params.pangenome_progressive
    (default false). Reuses the cactus container (minigraph) via the 'cactus_tools' label.

    Input : tuple(taxid, names_str, fastas)   (space-joined PanSN names + reference-first fastas)
    Output: growth = <taxid>.progressive_growth.tsv (k, sample_added, nodes, edges, bp) / versions
========================================================================================
*/

process PANGENOME_PROGRESSIVE {
    tag "taxid_${taxid}"
    label 'cactus_tools'

    publishDir "${params.outdir}/pangenome/${taxid}/progressive", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(names_str), path(fastas)

    output:
    tuple val(taxid), path("${taxid}.progressive_growth.tsv"), emit: growth
    path("versions.tsv"),                                      emit: versions

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"
    NAMES=(${names_str})
    FASTAS=(${fastas})
    N=\${#FASTAS[@]}

    # nodes, edges, bp of a minigraph rGFA (bp from S-line sequence, or LN:i: tag if seq is '*')
    measure() {
        awk 'BEGIN{FS="\\t"}
             \$1=="S"{ n++; if(\$3=="*"){ for(i=4;i<=NF;i++) if(\$i ~ /^LN:i:/){ v=\$i; sub(/^LN:i:/,"",v); bp+=v } } else bp+=length(\$3) }
             \$1=="L"{ e++ }
             END{ printf "%d\\t%d\\t%.0f\\n", n+0, e+0, bp+0 }' "\$1"
    }

    printf 'k\\tsample_added\\tnodes\\tedges\\tbp\\n' > ${taxid}.progressive_growth.tsv

    # k=1: reference-only graph
    minigraph -cxggs -t ${task.cpus} "\${FASTAS[0]}" > graph.gfa
    read n e bp < <(measure graph.gfa)
    printf '1\\t%s\\t%s\\t%s\\t%s\\n' "\${NAMES[0]}" "\$n" "\$e" "\$bp" >> ${taxid}.progressive_growth.tsv

    # add one assembly at a time, measuring after each
    for ((i=1; i<N; i++)); do
        minigraph -cxggs -t ${task.cpus} graph.gfa "\${FASTAS[\$i]}" > graph_new.gfa
        mv graph_new.gfa graph.gfa
        read n e bp < <(measure graph.gfa)
        printf '%d\\t%s\\t%s\\t%s\\t%s\\n' "\$((i+1))" "\${NAMES[\$i]}" "\$n" "\$e" "\$bp" >> ${taxid}.progressive_growth.tsv
    done
    rm -f graph.gfa

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tminigraph\\t%s\\n' "${task.process}" "\$(minigraph --version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    printf 'k\\tsample_added\\tnodes\\tedges\\tbp\\n1\\tref\\t1\\t0\\t0\\n' > ${taxid}.progressive_growth.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
