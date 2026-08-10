/*
========================================================================================
    ODGI STATS -> MultiQC MODULE  (reporting / workstream C)
========================================================================================
    Repo location: modules/pangenome_odgi_stats_mqc.nf

    Emits odgi graph statistics in the exact YAML format MultiQC's odgi module ingests
    (`odgi stats -i x.og -m -sgdl`), which MultiQC discovers via the "*.og.stats.yaml"
    filename pattern. One YAML for the whole-genome graph plus one per chromosome, so the
    MultiQC report shows a whole-genome row and a per-chromosome comparison table (nodes,
    edges, paths, per-component acyclicity, self-loops, nucleotide composition, mean links
    length, sum of path-node distances).

    NB the MultiQC odgi module requires EXACTLY the `-m -s -g -d -l` output (flags separate --
    this build rejects bundled short flags), and `odgi stats -m` requires OPTIMIZED node IDs, so
    each graph is `odgi sort -O`'d on a throwaway copy first.

    Input : tuple(taxid, og, [chrom_ogs])   (clip whole .og + per-chrom clip/full .og list)
    Output: yaml (whole + per-chrom *.og.stats.yaml) / versions
========================================================================================
*/

process PANGENOME_ODGI_STATS_MQC {
    tag "taxid_${taxid}"
    label 'cactus_tools'

    publishDir "${params.outdir}/pangenome/${taxid}/multiqc_stats", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(og), path(chrom_ogs)

    output:
    tuple val(taxid), path("*.og.stats.yaml"), emit: yaml
    path("versions.tsv"),                      emit: versions

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"

    # odgi stats -m requires OPTIMIZED (compacted) node IDs; cactus's .og is not, so sort -O
    # a throwaway copy of each graph first (the flags must also be separate -- this build
    # rejects bundled short flags like -sgdl).
    if odgi sort -i ${og} -o whole.sorted.og -O 2>/dev/null; then
        odgi stats -i whole.sorted.og -m -s -g -d -l > ${taxid}.og.stats.yaml || echo "whole-genome odgi stats failed" >&2
        rm -f whole.sorted.og
    else
        echo "whole-genome odgi sort failed" >&2
    fi

    # per-chromosome clip graphs -> one YAML each (sample name = chrN_1)
    for o in ${chrom_ogs}; do
        case "\$o" in *.full.og) continue;; esac
        base=\$(basename "\$o" .og)
        if odgi sort -i "\$o" -o "\${base}.sorted.og" -O 2>/dev/null; then
            odgi stats -i "\${base}.sorted.og" -m -s -g -d -l > "\${base}.og.stats.yaml" || echo "odgi stats failed: \$base" >&2
            rm -f "\${base}.sorted.og"
        else
            echo "odgi sort failed: \$base" >&2
        fi
    done

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\todgi\\t%s\\n' "${task.process}" "\$(odgi version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    printf -- '---\\nlength: 0\\nnodes: 0\\nedges: 0\\npaths: 0\\n' > ${taxid}.og.stats.yaml
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
