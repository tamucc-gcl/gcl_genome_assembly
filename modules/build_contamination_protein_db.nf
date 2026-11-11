process BUILD_CONTAMINATION_PROTEIN_DB {
    tag "${release}"
    label 'contamination_protein_db'
    conda 'bioconda::diamond=2.1.14 conda-forge::python=3.11 conda-forge::curl'
    cpus 16
    memory '64 GB'
    time '48h'
    // Build on shared work storage so the database does not rely on node-local stage-out.
    scratch false
    storeDir { "${db_root}/diamond/uniprot_reference_proteomes/${release}" }

    input:
    val db_root
    val release
    val source
    val source_sha256
    val source_format
    val taxonomy_dir
    path builder

    output:
    tuple path('reference_proteomes.dmnd'), path('protein_database.json'), emit: database

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    python3 ${quote(builder)} --source ${quote(source)} --source-sha256 ${quote(source_sha256)} \
      --format ${quote(source_format)} --release ${quote(release)} --taxonomy ${quote(taxonomy_dir)} \
      --threads ${task.cpus}
    """
}
