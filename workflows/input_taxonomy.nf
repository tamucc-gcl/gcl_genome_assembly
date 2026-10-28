process RESOLVE_SCIENTIFIC_NAMES {
    tag 'scientific_names'
    cpus 1
    memory '2 GB'
    time '1h'
    conda 'conda-forge::python=3.11'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/taxonomy", mode: params.publish_dir_mode
    input:
    val requests
    path names
    path nodes
    path resolver
    output:
    path 'resolved_names.tsv', emit: identities
    script:
    def payload=groovy.json.JsonOutput.toJson(requests)
    """
    set -euo pipefail
    cat > requests.json <<'NAME_REQUESTS'
${payload}
NAME_REQUESTS
    python3 ${resolver} --requests requests.json --names '${names}' --nodes '${nodes}'
    """
}

workflow INPUT_TAXONOMY {
    take:
    input_records
    taxdump
    main:
    input_records.map { meta, reads ->
        if (meta.taxid != null && !(meta.taxid.toString() ==~ /[1-9][0-9]*/))
            error "Invalid taxid for '${meta.sample}': expected a positive integer"
        if (meta.taxid == null && !meta.species?.trim())
            error "Missing organism identity for '${meta.sample}': supply taxid or species"
        tuple(meta, reads)
    }.branch { meta, reads ->
        known: meta.taxid != null
        names: true
    }.set { routes }
    requests=routes.names.toList().filter { !it.isEmpty() }
        .map { rows -> rows.collect { meta, reads -> [sample:meta.sample, species:meta.species] } }
    RESOLVE_SCIENTIFIC_NAMES(requests,
        taxdump.map { dir -> file("${dir}/names.dmp",checkIfExists:true) },
        taxdump.map { dir -> file("${dir}/nodes.dmp",checkIfExists:true) },
        file("${projectDir}/py_scripts/resolve_taxonomy_names.py",checkIfExists:true))
    resolved=RESOLVE_SCIENTIFIC_NAMES.out.identities.splitCsv(header:true,sep:'\t')
        .map { row -> tuple(row.sample.toString(),row.taxid.toString(),row.species.toString()) }
    named=routes.names.map { meta,reads -> tuple(meta.sample,meta,reads) }
        .join(resolved, failOnDuplicate:true, failOnMismatch:true)
        .map { sample,meta,reads,taxid,species -> tuple(meta+[taxid:taxid,species:species],reads) }
    emit:
    samples=routes.known.mix(named)
}
