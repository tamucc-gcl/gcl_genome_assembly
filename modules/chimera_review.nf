process CHIMERA_REVIEW {
    tag "${meta.id}"
    label 'chimera_review'
    publishDir "${params.outdir}/assembly/chimeras/review", mode: params.publish_dir_mode
    input:
    tuple val(meta), path(calls), path(context), path(audit), path(name_map), path(supplement,stageAs: 'supplement/*')
    path(script)
    output:
    tuple val(meta), path("${meta.id}.review/review.tsv"), emit: rows
    tuple val(meta), path("${meta.id}.review"), emit: packets
    script:
    """
    python3 chimera_review.py --calls '${calls}' --assembly '${meta.id}' --context '${context}' --coordinate-audit '${audit}' --name-map '${name_map}' --supplement supplement --supplement-status ${params.chimera_evidence == false ? 'disabled' : 'enabled'} --out '${meta.id}.review'
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
    path('README.md'), emit: report
    path('assembly_registry.tsv'), emit: registry
    path('cut-instructions.md'), emit: instructions
    script:
    def entryScript = script instanceof List ? script.find { it.name == 'chimera_review_index.py' } : script
    def args = (tables instanceof List ? tables : [tables]).collect { "--packet '${it}'" }.join(' ')
    """
    python3 ${entryScript} ${args}
    """
}
