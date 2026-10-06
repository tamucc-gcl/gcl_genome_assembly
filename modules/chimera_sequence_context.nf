process CHIMERA_SEQUENCE_CONTEXT {
    tag "${meta.id}"
    stageInMode 'symlink'
    label 'chimera_sequence_context'
    publishDir "${params.outdir}/assembly/chimeras/sequence_context", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta, stageAs: 'assessment/*'), path(calls),
          val(peer_meta), path(peer_fastas, stageAs: 'peer??/*'), path(reads, stageAs: 'reads/*'), val(motif),
          path(agp, stageAs: 'agp/*'), path(pairs, stageAs: 'pairs/*'),
          path(source, stageAs: 'source/*'), path(libraries, stageAs: 'libraries/*'),
          path(bams, stageAs:'reused/*'), path(bam_provenance, stageAs:'bam_provenance/*'), path(native_graph,stageAs:'graph/*')
    path(script)

    output:
    tuple val(meta), path("${meta.id}.sequence_context"), emit: context
    path('versions.tsv'), emit: versions

    script:
    def paths = peer_fastas instanceof List ? peer_fastas : [peer_fastas]
    def manifest = peer_meta.withIndex().collect { m, i -> [id: m.id, sample: m.sample, path: paths[i].toString(), auto_evidence:m.auto_evidence == true, chromosome_labels:m.chromosome_labels ?: [:]] }
    def peers = groovy.json.JsonOutput.toJson(manifest).replace("'", "'\"'\"'")
    def readArg = reads.toString().tokenize('/').last() == 'NO_PAIRS' ? '' : "--reads '${reads}'"
    def assayArg = agp.toString().tokenize('/').last() == 'NO_PAIRS' ? '' : "--agp '${agp}' --pairs '${pairs}' --pairs-source '${source}' --libraries '${libraries}'"
    def entryScript = script instanceof List ? script.find { it.name == 'chimera_sequence_context.py' } : script
    def bamPaths = bams instanceof List ? bams : [bams]
    def bam = bamPaths.find { it.name.endsWith('.bam') }
    def bamArg = bam ? "--bam '${bam}' --bam-provenance '${bam_provenance}'" : ''
    def graphArg = native_graph.toString().tokenize('/').last() == 'NO_PAIRS' ? '' : "--native-graph '${native_graph}'"
    """
    set -euo pipefail
    python3 ${entryScript} --fasta '${fasta}' --calls '${calls}' \\
        --assembly '${meta.id}' --sample '${meta.sample}' --peers '${peers}' \\
        --motif '${motif}' --threads ${task.cpus} ${readArg} ${bamArg} ${assayArg} ${graphArg} --out '${meta.id}.sequence_context'
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    for tool in minimap2 samtools tidk; do
        printf '${task.process}\\t%s\\t%s\\n' "\$tool" "\$(\$tool --version | sed -n '1p')" >> versions.tsv
    done
    """
}
