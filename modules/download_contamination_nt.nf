process DOWNLOAD_CONTAMINATION_NT {
    tag "${collection}/${release}"
    conda 'bioconda::blast=2.16.0 conda-forge::python=3.11 conda-forge::curl'
    scratch false
    cpus 2
    memory '4 GB'
    time '72h'
    storeDir { "${db_root}/blast/${collection}/${release}" }

    input:
    val db_root
    val collection
    val release
    path downloader

    output:
    path 'database', emit: bundle

    script:
    """
    python3 download_contamination_nt.py download --collection '${collection}' \
        --release '${release}' --out database
    """
}

process VERIFY_CONTAMINATION_NT {
    conda 'conda-forge::python=3.11'
    memory '1 GB'
    time '30m'
    cache false

    input:
    path bundle
    val db_root
    val collection
    val release
    path downloader

    output:
    val prefix, emit: prefix
    path 'nucleotide_database.json', emit: manifest

    script:
    prefix = "${db_root}/blast/${collection}/${release}/database/${collection}"
    """
    python3 download_contamination_nt.py verify --collection '${collection}' \
        --release '${release}' --out database
    cp database/database.json nucleotide_database.json
    """
}
