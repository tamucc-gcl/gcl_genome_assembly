nextflow.enable.dsl = 2

params.taxonomy_dir = null
params.contamination_protein_db_root = null
params.contamination_protein_release = null
params.contamination_protein_source = null
params.contamination_protein_sha256 = null
params.contamination_protein_format = 'tar'
params.contamination_protein_taxonomy_dir = null

include { CONTAMINATION_PROTEIN_DB } from '../../workflows/contamination_protein_db.nf'

workflow {
    if (!params.contamination_protein_db_root || !params.contamination_protein_release)
        error 'Supply --contamination_protein_db_root and --contamination_protein_release'
    def root = file(params.contamination_protein_db_root).toAbsolutePath().toString()
    def requestedTaxonomy = params.contamination_protein_taxonomy_dir ?: params.taxonomy_dir
    def taxonomy = requestedTaxonomy ? file(requestedTaxonomy).toAbsolutePath().toString() : null
    def source = params.contamination_protein_source
    if (source && !source.toString().startsWith('https://'))
        source = file(source).toAbsolutePath().toString()
    CONTAMINATION_PROTEIN_DB([root: root, release: params.contamination_protein_release.toString(),
        source: source, sha256: params.contamination_protein_sha256, format: params.contamination_protein_format,
        taxonomy: taxonomy],
        file("${projectDir}/../..").normalize())
    CONTAMINATION_PROTEIN_DB.out.directory.view { directory -> "Protein database bundle: ${directory}" }
}
