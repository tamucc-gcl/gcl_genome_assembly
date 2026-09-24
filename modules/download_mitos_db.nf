// Each requested RefSeq directory is a mandatory storeDir output.
process DOWNLOAD_MITOS_DB {
    tag "mitos_db"
    label 'mitos_download'
    storeDir { db_dir }

    input:
    val db_dir
    val refseq_sets
    val force

    output:
    val db_dir, emit: db
    path(refseq_sets.tokenize(',').collect { it.trim() }), emit: datasets

    script:
    def sets = refseq_sets.tokenize(',').collect { it.trim() }
    if (!sets || sets.any { !(it ==~ /refseq[0-9]+[mfo]/) })
        error "Invalid mitos_refseq_sets: ${refseq_sets}"
    def downloads = sets.collect { s ->
        """
        curl --fail --location --retry 3 -o ${s}.tar.bz2 \\
            'https://zenodo.org/records/3685310/files/${s}.tar.bz2?download=1'
        tar -xjf ${s}.tar.bz2
        test -d ${s}
        test -n "\$(find ${s} -type f -size +0c -print -quit)"
        """
    }.join('\n')
    """
    set -euo pipefail
    ${downloads}
    """

    stub:
    """
    echo 'Stub database downloads are disabled; use a prepared test database.' >&2
    exit 1
    """
}
