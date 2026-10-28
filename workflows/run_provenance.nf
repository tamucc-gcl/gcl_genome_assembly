include { COLLECT_SOFTWARE_VERSIONS } from '../modules/collect_software_versions.nf'

workflow RUN_PROVENANCE {
    take:
    versions
    main:
    sources=versions.toList().map { files -> files.unique { it.toString() }.sort { it.toString() } }
    def catalogs=[file("${projectDir}/nextflow.config"),file("${projectDir}/main.nf")]+
        file("${projectDir}/environments/*.{yaml,yml}")+file("${projectDir}/modules/*.nf")+file("${projectDir}/workflows/*.nf")+
        file("${projectDir}/functions/*.nf")+file("${projectDir}/py_scripts/*.py")+file("${projectDir}/r_scripts/*.R")
    catalogs=catalogs.sort { it.toString() }
    def context=[pipeline:workflow.manifest.name?.toString(),revision:workflow.revision?.toString(),
        commit:workflow.commitId?.toString(),run_name:workflow.runName.toString(),
        nextflow:nextflow.version.toString(),profile:workflow.profile.toString(),start:workflow.start.toString(),
        parameters:params.findAll { key,value -> !(key.toString() =~ /(?i)(password|secret|credential|token)/) },
        catalog_sources:catalogs.collect { it.toString() }]
    COLLECT_SOFTWARE_VERSIONS(sources,sources.map { files -> files.collect { [path:it.toString()] } },
        context,file("${projectDir}/py_scripts/run_provenance.py",checkIfExists:true),catalogs)
    emit:
    versions=COLLECT_SOFTWARE_VERSIONS.out.versions
    files=COLLECT_SOFTWARE_VERSIONS.out.files
}
