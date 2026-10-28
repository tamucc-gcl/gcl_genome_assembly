include { comparisonPlan } from '../functions/comparison_pairs.nf'
include { PAIRWISE_ALIGNMENT; PAIRWISE_PLOTS; COLLECT_PAIRWISE_RESULTS; COMPARISON_PLAN } from '../modules/pairwise_alignment.nf'

workflow PAIRWISE_COMPARISONS {
    take:
    final_assemblies // (meta, final FASTA, final FAI)
    main:
    helper=file("${projectDir}/py_scripts/pairwise_paf.py",checkIfExists:true)
    plot_helpers=[helper,file("${projectDir}/py_scripts/plot_pairwise.py",checkIfExists:true)]
    cohort=final_assemblies.toList()
    plan=cohort.map { rows -> comparisonPlan(rows,params.pairwise_alignment_mode ?: 'within_sample') }
    COMPARISON_PLAN(plan.map { it.inventory },cohort.map { rows -> rows.collect { it[2] } },
        cohort.map { rows -> rows.collect { m,fa,fai -> [id:m.id,sample:m.sample,taxid:m.taxid.toString(),
            fasta:fa.toString(),fai:fai.toString(),index_name:fai.name] } })
    pairs=plan.flatMap { it.pairs }
    PAIRWISE_ALIGNMENT(pairs,helper)
    plot_inputs=PAIRWISE_ALIGNMENT.out.paf.map { a,b,paf -> tuple([a,b],paf) }
        .join(pairs.map { a,fa1,fai1,b,fa2,fai2 -> tuple([a,b],fai1,fai2) },failOnDuplicate:true,failOnMismatch:true)
        .map { pair,paf,fai1,fai2 -> tuple(pair[0],pair[1],paf,fai1,fai2) }
    PAIRWISE_PLOTS(plot_inputs,plot_helpers)
    COLLECT_PAIRWISE_RESULTS(PAIRWISE_ALIGNMENT.out.qc.map { a,b,qc -> qc }.toList(),helper)
    emit:
    plan=COMPARISON_PLAN.out.plan
    inventory=COMPARISON_PLAN.out.inventory
    summary=COLLECT_PAIRWISE_RESULTS.out.summary
    dotplot=PAIRWISE_PLOTS.out.dotplot
    riparian=PAIRWISE_PLOTS.out.riparian
    status=PAIRWISE_PLOTS.out.status
    versions=PAIRWISE_ALIGNMENT.out.versions.mix(PAIRWISE_PLOTS.out.versions)
}
