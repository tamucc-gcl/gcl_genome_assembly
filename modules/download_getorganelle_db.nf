// Nextflow moves completed task-local database directories into the shared store.
process DOWNLOAD_GETORGANELLE_DB {
    tag "getorganelle_db"
    label 'getorganelle_download'
    storeDir { db_dir }

    input:
    val db_dir
    val force_download

    output:
    val db_dir, emit: config_dir
    path 'SeedDatabase', emit: seeds
    path 'LabelDatabase', emit: labels

    script:
    """
    set -euo pipefail
    get_organelle_config.py -a all --config-dir . --verbose
    # The downloader can continue after a failed individual database fetch.
    for kind in embplant_pt embplant_mt embplant_nr fungus_mt fungus_nr animal_mt other_pt; do
        test -s "SeedDatabase/\$kind.fasta"
        test -s "LabelDatabase/\$kind.fasta"
    done
    test -s SeedDatabase/VERSION
    test -s LabelDatabase/VERSION
    get_organelle_config.py --check --config-dir .
    """

    stub:
    """
    echo 'Stub database downloads are disabled; use a prepared test database.' >&2
    exit 1
    """
}
