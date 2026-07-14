/*
========================================================================================
    MITO CIRCULAR MAP MODULE
========================================================================================
    Generates a circular gene map of the annotated mitochondrial genome
    using pyCirclize, from the GenBank annotation produced by MITOHIFI.
    Repo location: modules/mito_circular_map.nf

    Input:
    - meta + genbank: Annotated mitogenome GenBank file from MITOHIFI
    Output:
    - circular_map: Publication-quality circular gene map PNG

    Dependencies: Python 3, pycirclize, biopython, matplotlib
========================================================================================
*/

process MITO_CIRCULAR_MAP {
    tag "${meta.sample}"
    label 'mito_circular_map'

    publishDir "${params.outdir}/mitogenome", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(genbank)
    path(mito_circular_script)

    output:
    tuple val(meta), path("${meta.sample}_mito_circular.png"), emit: circular_map

    script:
    def dpi = params.mitohifi_circular_map_dpi ?: 300
    """
    set -euo pipefail

    python3 ${mito_circular_script} \\
        --genbank ${genbank} \\
        --sample_id ${meta.sample} \\
        --output ${meta.sample}_mito_circular.png \\
        --dpi ${dpi}
    """

    stub:
    """
    touch ${meta.sample}_mito_circular.png
    """
}
