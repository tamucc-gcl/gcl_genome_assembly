/* Cohort-wide discovery, one HiFi map per source assembly, evidence, and review. */
process MISASSEMBLY_CATALOG {
    label 'misassembly_discovery'
    input:
    tuple val(records), path(name_maps,stageAs:'maps??/*'), val(pairs), path(pafs,stageAs:'alignments??/*')
    path(helpers)
    output:
    path('catalog.json'), emit: catalog
    script:
    def maps = name_maps instanceof List ? name_maps : [name_maps]
    def paths = pafs instanceof List ? pafs : [pafs]
    def manifest = records.withIndex().collect { m,i -> m + [name_map:maps[i].toString()] }
    def json = groovy.json.JsonOutput.toJson(manifest).replace("'", "'\"'\"'")
    def alignments = groovy.json.JsonOutput.toJson(pairs.withIndex().collect { p,i -> [a:p.a,b:p.b,path:paths[i].toString()] }).replace("'", "'\"'\"'")
    """
    python3 misassembly_discover.py catalog --manifest '${json}' --alignments '${alignments}' --out catalog.json
    """
}

process MISASSEMBLY_ALIGN {
    tag "${pair.a} vs ${pair.b}"
    label 'misassembly_alignment'
    input:
    tuple val(pair), path(query,stageAs:'query/*'), path(target,stageAs:'target/*')
    output:
    tuple val(pair), path("${pair.key}.paf.gz"), emit: alignment
    script:
    """
    set -euo pipefail
    samtools faidx '${target}'
    awk '{n+=\$2} END {if (n>=8000000000) exit 2}' '${target}.fai'
    minimap2 -x asm5 -I 8G -c --eqx --secondary=yes -N 20 -p 0.5 -t ${task.cpus} '${target}' '${query}' | gzip -c > '${pair.key}.paf.gz'
    """
}

process MISASSEMBLY_DISCOVER {
    tag "${meta.id}"
    label 'misassembly_discovery'
    input:
    tuple val(meta), path(fasta,stageAs:'assembly/*'), path(catalog), val(pairs), path(pafs,stageAs:'alignments??/*')
    path(helpers)
    output:
    tuple val(meta), path("${meta.id}.discovery"), emit: discovery
    script:
    def paths = pafs instanceof List ? pafs : [pafs]
    def manifest = pairs.withIndex().collect { p,i -> [a:p.a,b:p.b,path:paths[i].toString()] }
    def json = groovy.json.JsonOutput.toJson(manifest).replace("'", "'\"'\"'")
    """
    python3 misassembly_discover.py discover --catalog '${catalog}' --assembly '${meta.id}' --fasta '${fasta}' --alignments '${json}' --minimum-arm-bp ${params.chimera_min_arm_bp} --workers ${task.cpus} --out '${meta.id}.discovery'
    """
}

process MISASSEMBLY_MAP_HIFI {
    tag "${meta.id}"
    label 'misassembly_assessment'
    input:
    tuple val(meta), path(fasta,stageAs:'assembly/*'), path(discovery), path(reads,stageAs:'reads/*'), path(reused,stageAs:'reused/*'), path(provenance,stageAs:'provenance/*')
    path(helpers)
    output:
    tuple val(meta), path("${meta.id}.mapping"), emit: mapping
    script:
    def bamPaths = reused instanceof List ? reused : [reused]
    def bam = bamPaths.find { it.name.endsWith('.bam') }
    def index = bamPaths.find { it.name.endsWith('.bai') }
    def source = bam ? "--bam '${bam}' --index '${index}' --provenance '${provenance}'" : reads.name.endsWith('NO_PAIRS') ? '' : "--reads '${reads}'"
    """
    python3 misassembly_map.py --fasta '${fasta}' --discovery '${discovery}' ${source} --threads ${task.cpus} --out '${meta.id}.mapping'
    """
}

process MISASSEMBLY_ASSESS {
    tag "${meta.id}"
    label 'misassembly_assessment'
    input:
    tuple val(meta), path(fasta,stageAs:'assembly/*'), path(discovery), path(mapping), path(agp,stageAs:'agp/*'), path(pairs,stageAs:'pairs/*'), path(source,stageAs:'source/*'), path(libraries,stageAs:'libraries/*'), path(graph,stageAs:'graph/*'), val(motif)
    path(helpers)
    output:
    tuple val(meta), path("${meta.id}.evidence"), emit: evidence
    script:
    def has = { p -> !p.name.endsWith('NO_PAIRS') }
    def agpArg = has(agp) ? "--agp '${agp}'" : ''
    def hicArg = has(pairs) ? "--pairs '${pairs}' --source '${source}' --libraries '${libraries}'" : ''
    def graphArg = has(graph) ? "--graph '${graph}'" : ''
    """
    python3 misassembly_assess.py --fasta '${fasta}' --discovery '${discovery}' --mapping '${mapping}/mapping.json' ${agpArg} ${hicArg} ${graphArg} --motif '${motif}' --threads ${task.cpus} --plots ${params.chimera_plots} --out '${meta.id}.evidence'
    """
}

process MISASSEMBLY_REVIEW {
    label 'chimera_joins'
    publishDir "${params.outdir}/assembly/chimeras", mode: params.publish_dir_mode
    input:
    path(packets,stageAs:'packets??/*')
    path(helpers)
    path(actions,stageAs:'actions/*')
    path(verifications,stageAs:'verified??/*')
    output:
    path('README.md'), emit: report
    path('review.tsv'), emit: decisions
    path('assembly_registry.tsv'), emit: registry
    path('mapping_manifest.json'), emit: mapping_manifest
    path('cut-instructions.md'), emit: instructions
    path('coordinate_status.md'), emit: coordinates
    path('assemblies'), emit: assemblies
    path('versions.tsv'), emit: versions
    script:
    def args = (packets instanceof List ? packets : [packets]).collect { "--packet '${it}'" }.join(' ')
    def actionArg = actions.name.endsWith('NO_PAIRS') ? '' : "--actions '${actions}'"
    def verificationArgs = (verifications instanceof List ? verifications : [verifications]).findAll { !it.name.endsWith('NO_PAIRS') }.collect { "--verification '${it}'" }.join(' ')
    """
    python3 misassembly_report.py ${args} ${actionArg} ${verificationArgs}
    """
}
