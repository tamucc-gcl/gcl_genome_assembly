nextflow.enable.dsl = 2

params.contamination_db_root = null
params.proteome_source = null
params.proteome_sha256 = null
params.proteome_release = null
params.proteome_format = null
params.taxonomy_dir = null

include { BUILD_CONTAMINATION_PROTEIN_DB } from '../../modules/build_contamination_protein_db.nf'
include { DOWNLOAD_TAXDUMP } from '../../modules/download_taxdump.nf'

workflow {
    if (!params.contamination_db_root || !params.proteome_source ||
        !(params.proteome_sha256?.toString() ==~ /[a-fA-F0-9]{64}/) ||
        !(params.proteome_release?.toString() ==~ /[a-zA-Z0-9][a-zA-Z0-9_.-]*/) ||
        !(params.proteome_format in ['fasta', 'tar']))
        error 'Supply --contamination_db_root, --proteome_source, --proteome_sha256, --proteome_release and --proteome_format (fasta or tar)'

    def root = file(params.contamination_db_root).toAbsolutePath().toString()
    def source = params.proteome_source.toString()
    if (!source.startsWith('https://'))
        source = file(source, checkIfExists: true).toAbsolutePath().toString()
    def release = params.proteome_release.toString()
    def checksum = params.proteome_sha256.toString().toLowerCase()
    def format = params.proteome_format.toString()
    def taxonomy = params.taxonomy_dir ? file(params.taxonomy_dir).toAbsolutePath().toString() :
                                        "${root}/ncbi_taxonomy"
    def builder = file("${projectDir}/../../py_scripts/build_contamination_protein_db.py", checkIfExists: true)
    DOWNLOAD_TAXDUMP(taxonomy, false)
    BUILD_CONTAMINATION_PROTEIN_DB(root, release, source, checksum, format,
        DOWNLOAD_TAXDUMP.out.taxdump_dir, builder)
    BUILD_CONTAMINATION_PROTEIN_DB.out.database.map { db, manifest ->
        def identity = new groovy.json.JsonSlurper().parseText(manifest.text)
        if (identity.release != release || identity.source_sha256 != checksum ||
            identity.source_format != format || db.size() == 0)
            error "Stored protein database differs from the requested source. Choose a new release directory or repair ${manifest}."
        "Protein database: ${db}\nProvenance: ${manifest}"
    }.view()
}
