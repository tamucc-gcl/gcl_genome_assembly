include { DOWNLOAD_CONTAMINATION_PROTEIN_TAXONOMY; RESOLVE_CONTAMINATION_PROTEIN_SOURCE;
          DOWNLOAD_CONTAMINATION_PROTEIN_SOURCE } from '../modules/contamination_protein_source.nf'
include { BUILD_CONTAMINATION_PROTEIN_DB } from '../modules/build_contamination_protein_db.nf'

workflow CONTAMINATION_PROTEIN_DB {
    take:
    specification
    base

    main:
    def resolver = file(base.resolve('py_scripts/contamination_protein_source.py'), checkIfExists: true)
    def builder = file(base.resolve('py_scripts/build_contamination_protein_db.py'), checkIfExists: true)
    def taxonomyDirectory = specification.taxonomy ?: "${specification.root}/ncbi_taxonomy/${specification.release}"
    DOWNLOAD_CONTAMINATION_PROTEIN_TAXONOMY(taxonomyDirectory.toString())
    def taxonomy = DOWNLOAD_CONTAMINATION_PROTEIN_TAXONOMY.out.directory
    RESOLVE_CONTAMINATION_PROTEIN_SOURCE(specification, resolver)
    def requests = RESOLVE_CONTAMINATION_PROTEIN_SOURCE.out.map { path ->
        tuple(path, new groovy.json.JsonSlurper().parseText(path.text))
    }.branch { path, request ->
        reuse: request.existing
        build: true
    }
    DOWNLOAD_CONTAMINATION_PROTEIN_SOURCE(specification.root, specification.release,
        requests.build.map { path, request -> path }, resolver)
    def downloaded = DOWNLOAD_CONTAMINATION_PROTEIN_SOURCE.out
        .combine(requests.build.map { path, request -> request })
        .map { archive, manifest, request ->
            def identity = new groovy.json.JsonSlurper().parseText(manifest.text)
            ['release', 'source', 'source_format'].each { field ->
                if (identity[field] != request[field]) error "Stored protein download differs from requested ${field}"
            }
            ['sha256', 'md5'].each { field ->
                if (request[field] && identity[field] != request[field])
                    error "Stored protein download checksum differs from requested ${field}"
            }
            if (archive.size() != identity.actual_bytes) error 'Stored protein archive size changed'
            tuple([location: archive.toString(), origin: request.source], identity.sha256, request.source_format)
        }
    BUILD_CONTAMINATION_PROTEIN_DB(specification.root, specification.release,
        downloaded.map { source, sha256, format -> source },
        downloaded.map { source, sha256, format -> sha256 },
        downloaded.map { source, sha256, format -> format }, taxonomy, builder)
    def existing = requests.reuse.map { path, request -> request.bundle.toString() }
    def built = BUILD_CONTAMINATION_PROTEIN_DB.out.database.map { db, manifest -> db.parent.toString() }
    def bundle = existing.mix(built).first()

    emit:
    directory = bundle
    taxdump = taxonomy
}
