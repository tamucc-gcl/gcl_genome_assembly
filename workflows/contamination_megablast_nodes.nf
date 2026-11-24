include { CONTAMINATION_MEGABLAST_PLAN; CONTAMINATION_MEGABLAST_REMAINING; CONTAMINATION_MEGABLAST_NODE;
          CONTAMINATION_MEGABLAST_COLLECT; CONTAMINATION_MEGABLAST_CHECK } from '../modules/contamination_nodes.nf'

workflow CONTAMINATION_MEGABLAST_NODES {
    take:
    prepared
    database
    base

    main:
    def helper = file(base.resolve('py_scripts/contamination_workers.py'), checkIfExists: true)
    def allQueries = prepared.map { id, metrics, queries -> tuple(id, queries) }.toList()
    CONTAMINATION_MEGABLAST_PLAN(
        allQueries.map { rows -> rows.collect { id, queries -> [assembly: id, queries: queries.toString()] } },
        allQueries.map { rows -> rows.collect { it[1] } }, database,
        workflow.workDir.resolve('contamination_megablast').toString(), helper)
    CONTAMINATION_MEGABLAST_REMAINING(CONTAMINATION_MEGABLAST_PLAN.out, helper)
    def requestedWorkers = CONTAMINATION_MEGABLAST_REMAINING.out.flatMap { status ->
        def remaining = new groovy.json.JsonSlurper().parseText(status.text)
        remaining.remaining ? (1..(params.contamination_megablast_nodes as int)).toList() : []
    }
    CONTAMINATION_MEGABLAST_NODE(requestedWorkers,
        CONTAMINATION_MEGABLAST_PLAN.out, helper)
    CONTAMINATION_MEGABLAST_COLLECT(prepared.map { id, metrics, queries -> id },
        CONTAMINATION_MEGABLAST_PLAN.out, helper)
    CONTAMINATION_MEGABLAST_CHECK(CONTAMINATION_MEGABLAST_PLAN.out,
        CONTAMINATION_MEGABLAST_NODE.out.toList(), helper)
    def hits = CONTAMINATION_MEGABLAST_COLLECT.out.flatMap { id, directory ->
        directory.toFile().listFiles().findAll { it.name.endsWith('.hits.tsv') }
            .sort { it.name }.collect { tuple(id, it.toPath()) }
    }

    emit:
    results = hits
}
