process CHIMERA_SEQUENCE_CONTEXT {
    tag "${meta.id}"
    label 'chimera_sequence_context'
    publishDir "${params.outdir}/assembly/chimeras/sequence_context", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta, stageAs: 'assessment/*'), path(calls),
          val(peer_meta), path(peer_fastas, stageAs: 'peer??/*'), path(reads, stageAs: 'reads/*'), val(motif)
    path(script)

    output:
    tuple val(meta), path("${meta.id}.sequence_context"), emit: context
    path('versions.tsv'), emit: versions

    script:
    def paths = peer_fastas instanceof List ? peer_fastas : [peer_fastas]
    def manifest = peer_meta.withIndex().collect { m, i -> [id: m.id, sample: m.sample, path: paths[i].toString()] }
    def peers = groovy.json.JsonOutput.toJson(manifest).replace("'", "'\"'\"'")
    def readArg = reads.toString().tokenize('/').last() == 'NO_PAIRS' ? '' : "--reads '${reads}'"
    """
    set -euo pipefail
    python3 ${script} --fasta '${fasta}' --calls '${calls}' \\
        --assembly '${meta.id}' --sample '${meta.sample}' --peers '${peers}' \\
        --motif '${motif}' --threads ${task.cpus} ${readArg} --out '${meta.id}.sequence_context'
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    for tool in minimap2 samtools tidk; do
        printf '${task.process}\\t%s\\t%s\\n' "\$tool" "\$(\$tool --version | sed -n '1p')" >> versions.tsv
    done
    """
}
