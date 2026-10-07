process CHIMERA_REVIEW {
    tag "${meta.id}"
    label 'chimera_joins'
    publishDir "${params.outdir}/assembly/chimeras/review", mode: params.publish_dir_mode
    input:
    tuple val(meta), path(calls), path(context)
    path(script)
    output:
    tuple val(meta), path("${meta.id}.review/review.tsv"), emit: rows
    tuple val(meta), path("${meta.id}.review"), emit: packets
    script:
    """
    python3 ${script} --calls '${calls}' --assembly '${meta.id}' --context '${context}' --out '${meta.id}.review'
    """
}

process CHIMERA_REVIEW_INDEX {
    label 'chimera_joins'
    publishDir "${params.outdir}/assembly/chimeras/review", mode: params.publish_dir_mode
    input:
    path(tables,stageAs:'review??/*')
    path(script)
    output:
    path('chimera_review.tsv'), emit: review_file
    path('index.html'), emit: report
    script:
    def args = (tables instanceof List ? tables : [tables]).collect { "--table '${it}'" }.join(' ')
    """
    python3 ${script} ${args}
    """
}
