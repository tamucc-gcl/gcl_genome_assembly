include { ASSEMBLY_SV_PLAN; ASSEMBLY_SV_ALIGN; ASSEMBLY_SV_CALL } from '../modules/assembly_sv.nf'

workflow ASSEMBLY_SV {
    take:
    eligibility
    main:
    reports = Channel.empty()
    versions = Channel.empty()
    if (params.run_assembly_sv) {
        ASSEMBLY_SV_PLAN(eligibility, file("${projectDir}/py_scripts/sv_plan.py", checkIfExists: true))
        reports = ASSEMBLY_SV_PLAN.out.report
        if (!params.assembly_sv_plan_only) {
            pairs = ASSEMBLY_SV_PLAN.out.plan.flatMap { plan ->
                new groovy.json.JsonSlurper().parseText(plan.text).pairs.collect { p ->
                    tuple(p, file(p.reference_fasta, checkIfExists: true), file(p.query_fasta, checkIfExists: true))
                }
            }
            ASSEMBLY_SV_ALIGN(pairs)
            ASSEMBLY_SV_CALL(ASSEMBLY_SV_ALIGN.out.alignment,
                file("${projectDir}/py_scripts/sv_alignment_qc.py", checkIfExists: true))
            reports = reports.mix(ASSEMBLY_SV_CALL.out.report.map { meta, report -> report })
            versions = ASSEMBLY_SV_ALIGN.out.versions.mix(ASSEMBLY_SV_CALL.out.versions)
        }
    }
    emit:
    report = reports
    versions = versions
}
