def nodeQuote(value) {
    "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process CONTAMINATION_MEGABLAST_PLAN {
    conda 'conda-forge::python=3.11 bioconda::blast=2.16.0'
    cpus 1
    memory '2 GB'
    time '2h'

    input:
    val specifications
    path sources, stageAs: 'sources/source??/*'
    val database
    val state
    path helper

    output:
    path 'plan.json'

    script:
    """
    cat > specifications.json <<'SPECIFICATIONS'
${groovy.json.JsonOutput.toJson(specifications)}
SPECIFICATIONS
    python3 contamination_workers.py plan --specs specifications.json \
        --database ${nodeQuote(database)} --state ${nodeQuote(state)} \
        --threads ${params.contamination_megablast_threads} --out plan.json
    """
}

process CONTAMINATION_MEGABLAST_REMAINING {
    conda 'conda-forge::python=3.11'
    cpus 1
    memory '1 GB'
    time '1h'
    cache false

    input:
    path plan
    path helper

    output:
    path 'remaining.json'

    script:
    """
    python3 contamination_workers.py remaining --plan ${nodeQuote(plan)} --out remaining.json
    """
}

process CONTAMINATION_MEGABLAST_POOL {
    // The coordinator owns its Slurm children so unused pending requests can retire.
    executor 'local'
    conda 'conda-forge::python=3.11 bioconda::blast=2.16.0'
    cpus 1
    memory '1 GB'
    scratch false
    cache false
    errorStrategy 'terminate'
    time { params.contamination_megablast_time }
    publishDir "${params.outdir}/pipeline/contamination_workers", mode: 'copy'
    input:
    path plan
    path helper
    path pool_helper
    output:
    path 'pool', emit: pool
    path 'pool/worker_*.json', optional: true, emit: summaries
    script:
    def seconds = nextflow.util.Duration.of(params.contamination_megablast_time.toString()).toSeconds()
    def walltime = "${(seconds / 3600) as int}:${((seconds % 3600) / 60) as int}:${seconds % 60}"
    """
    python3 ${nodeQuote(pool_helper)} --plan ${nodeQuote(plan)} --helper ${nodeQuote(helper)} --out pool \
        --nodes ${params.contamination_megablast_nodes} --slots ${params.contamination_megablast_slots} \
        --partition ${nodeQuote(params.contamination_megablast_partition)} \
        --nodelist ${nodeQuote(params.contamination_megablast_nodelist)} \
        --scratch ${nodeQuote(params.contamination_megablast_scratch)} --walltime ${nodeQuote(walltime)} \
        --timeout ${seconds}
    """
}

process CONTAMINATION_MEGABLAST_NODE {
    tag "node_worker_${worker}"
    conda 'conda-forge::python=3.11 bioconda::blast=2.16.0'
    cpus { (params.contamination_megablast_slots as int) * (params.contamination_megablast_threads as int) }
    // Nextflow 23.10 omits its memory flag for zero; request all node RAM explicitly.
    memory '0 B'
    clusterOptions '--exclusive --mem=0'
    scratch false
    maxForks (params.contamination_megablast_nodes as int)
    time { params.contamination_megablast_time }

    input:
    val worker
    path plan
    path helper

    output:
    path 'worker_summary.json'

    script:
    """
    python3 contamination_workers.py worker --plan ${nodeQuote(plan)} --worker ${worker} \
        --slots ${params.contamination_megablast_slots} \
        --scratch-root ${nodeQuote(params.contamination_megablast_scratch)} --out worker_summary.json
    """
}

process CONTAMINATION_MEGABLAST_COLLECT {
    tag "${assembly}"
    conda 'conda-forge::python=3.11'
    cpus 1
    memory '1 GB'
    time { params.contamination_megablast_time }

    input:
    val assembly
    path plan
    path helper

    output:
    tuple val(assembly), path('hits')

    script:
    def timeout = nextflow.util.Duration.of(params.contamination_megablast_time.toString()).toSeconds()
    """
    python3 contamination_workers.py collect --plan ${nodeQuote(plan)} --assembly ${nodeQuote(assembly)} \
        --timeout ${timeout} --out hits
    """
}

process CONTAMINATION_MEGABLAST_CHECK {
    conda 'conda-forge::python=3.11'
    cpus 1
    memory '1 GB'
    time '1h'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/pipeline", mode: 'copy', pattern: 'chunk_completion.json',
        saveAs: { 'contamination_megablast_completion.json' }

    input:
    path plan
    path summaries, stageAs: 'workers/worker??/*'
    path helper

    output:
    path 'chunk_completion.json'

    script:
    def files = summaries instanceof List ? summaries.collect { nodeQuote(it) }.join(' ') : nodeQuote(summaries)
    """
    python3 contamination_workers.py check --plan ${nodeQuote(plan)} \
        --summaries ${files} --out chunk_completion.json
    """
}
