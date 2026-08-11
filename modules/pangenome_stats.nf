/*
========================================================================================
    PANGENOME STATS MODULE
========================================================================================
    Repo location: modules/pangenome_stats.nf

    Graph metrics from the cactus outputs (same container -- vg and odgi are bundled).
    vg stats on the GBZ (sequence/node/edge counts, path count); odgi stats on the .og.
    Both the clipped and full outputs may be present; we run on each file found.

    Input : tuple(taxid, gbz(list), og(list))
    Output: per-taxid stats text
========================================================================================
*/

process PANGENOME_STATS {
    tag "taxid_${taxid}"
    label 'cactus_tools'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(gbz), path(og)

    output:
    tuple val(taxid), path("${taxid}.vg_stats.txt"),   emit: vg_stats
    tuple val(taxid), path("${taxid}.odgi_stats.txt"), emit: odgi_stats, optional: true

    script:
    """
    set -euo pipefail

    : > ${taxid}.vg_stats.txt
    for g in ${gbz}; do
        echo "## \${g}" >> ${taxid}.vg_stats.txt
        vg stats -lz "\${g}" >> ${taxid}.vg_stats.txt
        echo >> ${taxid}.vg_stats.txt
    done

    if ls ${og} >/dev/null 2>&1; then
        : > ${taxid}.odgi_stats.txt
        for o in ${og}; do
            echo "## \${o}" >> ${taxid}.odgi_stats.txt
            odgi stats -i "\${o}" -S >> ${taxid}.odgi_stats.txt
            echo >> ${taxid}.odgi_stats.txt
        done
    fi
    """

    stub:
    """
    : > ${taxid}.vg_stats.txt
    : > ${taxid}.odgi_stats.txt
    """
}
