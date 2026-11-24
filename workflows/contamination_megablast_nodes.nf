include { CONTAMINATION_MEGABLAST_PLAN; CONTAMINATION_MEGABLAST_REMAINING; CONTAMINATION_MEGABLAST_POOL;
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
    CONTAMINATION_MEGABLAST_POOL(CONTAMINATION_MEGABLAST_PLAN.out, helper,
        file(base.resolve('py_scripts/contamination_slurm_pool.py'),checkIfExists:true))
    // Collection follows pool completion; collectors no longer hold nodes while waiting.
    // combine a completion token with the plan rather than passing a channel inside a value.
    def readyPlan = CONTAMINATION_MEGABLAST_PLAN.out.combine(CONTAMINATION_MEGABLAST_POOL.out.pool).map { p,d -> p }
    CONTAMINATION_MEGABLAST_COLLECT(prepared.map { id, metrics, queries -> id },
        readyPlan.first(), helper)
    CONTAMINATION_MEGABLAST_CHECK(CONTAMINATION_MEGABLAST_PLAN.out,
        CONTAMINATION_MEGABLAST_POOL.out.summaries.flatMap { files -> files instanceof List ? files : [files] }.toList(), helper)
    def hits = CONTAMINATION_MEGABLAST_COLLECT.out.flatMap { id, directory ->
        directory.toFile().listFiles().findAll { it.name.endsWith('.hits.tsv') }
            .sort { it.name }.collect { tuple(id, it.toPath()) }
    }

    emit:
    results = hits
}
