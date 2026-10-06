nextflow.enable.dsl=2
include { CHIMERA } from './workflows/chimera.nf'

process VERIFY_CHIMERA_REAL {
    tag "${meta.id}"
    cpus 1
    memory '12 GB'
    publishDir "${params.outdir}/verification", mode:'copy'
    input:
    tuple val(meta), path(original,stageAs:'original/*'), path(corrected,stageAs:'corrected/*'), path(names)
    path(actions)
    path(scripts)
    output:
    path("${meta.id}.verification.json")
    script:
    def verifier = scripts.find { it.name == 'verify_chimera_real_test.py' }
    """
    python3 ${verifier} --original '${original}' --corrected '${corrected}' --names '${names}' \\
        --assembly '${meta.id}' --scenario '${params.scenario}' --actions '${actions}' --out '${meta.id}.verification.json'
    """
}

workflow {
    if (!params.manifest || !params.reviewed_actions_json) error 'Supply prepared manifest and reviewed action JSON'
    def data = new groovy.json.JsonSlurper().parseText(file(params.manifest,checkIfExists:true).text)
    def cohort = data.cohort
    def targets = data.targets
    CHIMERA(
        Channel.fromList(cohort.collect { r -> tuple(r.meta,file(r.fasta),file(r.names)) }),
        Channel.empty(),
        Channel.fromList(targets.collect { r -> tuple(r.meta,file(r.agp)) }),
        Channel.fromList(targets.collect { r -> tuple(r.meta.id,file(r.paf)) }),
        Channel.of(tuple(data.species,file(data.candidates))),
        Channel.of(tuple(data.species,file(data.reference_names))),
        Channel.of(tuple(data.species,data.reference_id)),
        Channel.of(file(data.quality)),
        Channel.fromList(targets.collect { r -> tuple(r.meta,'scaffold',file(r.pairs)) }),
        Channel.fromList(targets.collect { r -> tuple(r.meta,file(r.source),file(r.libraries)) }),
        Channel.fromList(targets.groupBy { it.meta.sample }.collect { sample,records -> tuple([sample:sample],records.collect { file(it.graph) }) }),
        Channel.of(tuple(data.species,'CCCTAA')),
        Channel.empty(),Channel.empty(),Channel.empty(),[scaffold:true,harmonize:true])
    original = Channel.fromList(cohort.collect { r -> tuple(r.meta.id,r.meta,file(r.fasta)) })
    VERIFY_CHIMERA_REAL(original.join(CHIMERA.out.pre_finalize.map { m,fa,nm -> tuple(m.id,fa,nm) })
        .map { id,m,source,corrected,names -> tuple(m,source,corrected,names) },
        file(params.reviewed_actions_json),
        Channel.value([file("${projectDir}/tests/verify_chimera_real_test.py"),
                       file("${projectDir}/py_scripts/break_chimeras.py"),file("${projectDir}/py_scripts/chimera_actions.py")]))
}
