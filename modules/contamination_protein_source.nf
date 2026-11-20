process DOWNLOAD_CONTAMINATION_PROTEIN_TAXONOMY {
    conda 'conda-forge::curl'
    cpus 1
    memory '1 GB'
    time '2h'
    scratch false
    storeDir { directory }
    input:
    val directory
    output:
    val directory, emit: directory
    path 'nodes.dmp'
    path 'names.dmp'
    path 'merged.dmp'
    script:
    """
    curl --fail --location --retry 3 -o taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzf taxdump.tar.gz nodes.dmp names.dmp merged.dmp
    test -s nodes.dmp
    test -s names.dmp
    test -s merged.dmp
    """
}

process RESOLVE_CONTAMINATION_PROTEIN_SOURCE {
    conda 'conda-forge::python=3.11 conda-forge::curl'
    cpus 1
    memory '1 GB'
    time '30m'
    cache false
    scratch false
    publishDir "${params.outdir}/pipeline", mode: 'copy', pattern: 'protein_request.json'

    input:
    val specification
    path resolver

    output:
    path 'protein_request.json'

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def custom = specification.source ? "--source ${quote(specification.source)} --sha256 ${quote(specification.sha256)}" : ''
    """
    python3 ${quote(resolver)} resolve --root ${quote(specification.root)} \
      --release ${quote(specification.release)} --format ${quote(specification.format)} \
      ${custom} --out protein_request.json
    """
}

process DOWNLOAD_CONTAMINATION_PROTEIN_SOURCE {
    tag "${release}"
    conda 'conda-forge::python=3.11 conda-forge::curl'
    cpus 2
    memory '2 GB'
    time '96h'
    scratch false
    storeDir { "${root}/downloads/uniprot_reference_proteomes/${release}" }

    input:
    val root
    val release
    path request
    path resolver

    output:
    tuple path('source/archive'), path('source/source.json')

    script:
    """
    python3 ${resolver} fetch --request '${request}' --out source
    """
}
