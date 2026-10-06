process CHIMERA_ADJUDICATE {
    tag "${meta.id}"
    label 'chimera_joins'
    publishDir "${params.outdir}/assembly/chimeras/decisions", mode: params.publish_dir_mode
    input:
    tuple val(meta), path(calls), path(context)
    path(script)
    output:
    tuple val(meta), path("${meta.id}.decisions/actions.tsv"), emit: actions
    tuple val(meta), path("${meta.id}.decisions"), emit: packets
    script:
    """
    python3 ${script} --calls '${calls}' --assembly '${meta.id}' --context '${context}' --out '${meta.id}.decisions' --min-piece-bp ${params.chimera_min_piece_bp ?: 1000000}
    """
}
