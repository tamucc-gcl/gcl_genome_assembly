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

    NB the MultiQC odgi module requires EXACTLY the `-m -sgdl` output; other flag
    combinations break it. Stats run on the graph as stored (no sort) -- topology metrics
    are order-independent and the linearity metrics reflect the graph's own node order.

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

    # whole-genome graph (sample name = species/taxid)
    odgi stats -i ${og} -m -sgdl > ${taxid}.og.stats.yaml || echo "whole-genome odgi stats failed" >&2

    # per-chromosome clip graphs -> one YAML each (sample name = chrN_1)
    for o in ${chrom_ogs}; do
        case "\$o" in *.full.og) continue;; esac
        base=\$(basename "\$o" .og)
        odgi stats -i "\$o" -m -sgdl > "\${base}.og.stats.yaml" || echo "odgi stats failed: \$base" >&2
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
