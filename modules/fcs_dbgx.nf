// Payload names documented in the NCBI FCS-GX quickstart.
def gxDatabaseFiles(manifest) {
    def prefix = manifest.replaceAll('.*/([^/]+)\\.manifest.*$', '$1')
    ['README.txt', 'assemblies.tsv', 'blast_div.tsv.gz', 'gxi', 'gxs',
     'manifest', 'meta.jsonl', 'seq_info.tsv.gz', 'taxa.tsv'].collect { "${prefix}.${it}" }
}

process FCS_DB_GET {
    tag "fcs_db_get"
    label 'fcs_download'
    storeDir { gxdb_dir }
    scratch false

    input:
    val gxdb_manifest
    val gxdb_dir
    val force_download

    output:
    val gxdb_dir, emit: out_dir
    path(gxDatabaseFiles(gxdb_manifest)), emit: files

    script:
    def required = gxDatabaseFiles(gxdb_manifest).collect { "test -s '${it}'" }.join('\n')
    """
    set -euo pipefail
    /app/bin/sync_files get --mft "${gxdb_manifest}" --dir .
    ${required}
    """

    stub:
    """
    echo 'Stub database downloads are disabled; use a prepared test database.' >&2
    exit 1
    """
}
