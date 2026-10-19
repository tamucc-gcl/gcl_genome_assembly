nextflow.enable.dsl=2

params.manifest = null
params.inputs = null
params.repo = null
params.outdir = null

process JOIN_LABELS {
    tag "${assembly}"
    cpus 1
    memory '32 GB'
    time '12h'
    input:
    tuple val(assembly), path(job)
    path helpers
    output:
    tuple val(assembly), path(job), path("${assembly}.labels"), emit: evidence
    script:
    """
    python3 '${params.repo}/scripts/comparisons/join_investigation.py' branch --kind labels \
      --job '${job}' --out '${assembly}.labels' --threads ${task.cpus}
    """
}

process JOIN_HIFI {
    tag "${assembly}"
    cpus 8
    memory '64 GB'
    time '24h'
    input:
    tuple val(assembly), path(job), path(labels)
    path helpers
    output:
    path("${assembly}.hifi"), emit: evidence
    script:
    """
    python3 '${params.repo}/scripts/comparisons/join_investigation.py' branch --kind hifi \
      --job '${job}' --masks '${labels}/masks.json' --out '${assembly}.hifi' --threads ${task.cpus}
    """
}

process JOIN_HIC {
    tag "${assembly}"
    cpus 1
    memory '32 GB'
    time '24h'
    input:
    tuple val(assembly), path(job), path(labels)
    path helpers
    output:
    path("${assembly}.hic"), emit: evidence
    script:
    """
    python3 '${params.repo}/scripts/comparisons/join_investigation.py' branch --kind hic \
      --job '${job}' --masks '${labels}/masks.json' --out '${assembly}.hic' --threads ${task.cpus}
    """
}

process JOIN_REPORT {
    cpus 1
    memory '12 GB'
    time '6h'
    publishDir params.outdir, mode: 'symlink'
    input:
    path packets
    path manifest
    path helpers
    output:
    path 'review'
    script:
    def args = packets.collect { "--packet '${it}'" }.join(' ')
    """
    python3 '${params.repo}/scripts/comparisons/join_investigation.py' report \
      --manifest '${manifest}' ${args} --out review
    """
}

workflow {
    if (!params.manifest || !params.inputs || !params.repo || !params.outdir)
        error 'Supply --manifest, --inputs, --repo and --outdir'
    def manifest = new groovy.json.JsonSlurper().parseText(file(params.manifest,checkIfExists:true).text)
    def helpers = Channel.value([
        file("${params.repo}/scripts/comparisons/join_investigation.py",checkIfExists:true),
        *['chimera_controls.py','chimera_tracks.py','chimera_intervals.py','misassembly_core.py'].collect {
            file("${params.repo}/py_scripts/${it}",checkIfExists:true)
        }
    ])
    def jobs = Channel.fromList(manifest.jobs).map { j ->
        tuple(j.assembly,file("${params.inputs}/${j.assembly}.json",checkIfExists:true))
    }
    JOIN_LABELS(jobs,helpers)
    JOIN_HIFI(JOIN_LABELS.out.evidence,helpers)
    JOIN_HIC(JOIN_LABELS.out.evidence,helpers)
    def packets = JOIN_LABELS.out.evidence.map { assembly,job,packet -> packet }
        .mix(JOIN_HIFI.out.evidence,JOIN_HIC.out.evidence).collect()
    JOIN_REPORT(packets,file(params.manifest,checkIfExists:true),helpers)
}
