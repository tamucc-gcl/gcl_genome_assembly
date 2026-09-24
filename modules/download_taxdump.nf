// Shared taxonomy snapshot: required files, rather than a sentinel, control reuse.
process DOWNLOAD_TAXDUMP {
    tag "download_taxdump"
    label 'download_db'
    storeDir { taxdump_dir }

    input:
    val taxdump_dir
    val force_download

    output:
    val taxdump_dir, emit: taxdump_dir
    path 'names.dmp', emit: names
    path 'nodes.dmp', emit: nodes

    script:
    """
    set -euo pipefail
    curl --fail --location --retry 3 -o taxdump.tar.gz \\
        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzf taxdump.tar.gz names.dmp nodes.dmp
    test -s names.dmp
    test -s nodes.dmp
    """

    stub:
    """
    echo 'Stub database downloads are disabled; use a prepared test database.' >&2
    exit 1
    """
}
