/*
========================================================================================
    DOWNLOAD GETORGANELLE DATABASES MODULE
========================================================================================
    Repo location: modules/download_getorganelle_db.nf

    Fetches ALL GetOrganelle seed + label databases once into a permanent, shared location
    before any parallel GETORGANELLE runs. Mirrors DOWNLOAD_TAXDUMP / FCS_DB_GET: download
    into an absolute dir, guard with a sentinel, skip if present, emit the config-dir path
    for downstream --config-dir. DBs are small (a few hundred MB total), so we grab
    everything (`-a all`). GetOrganelle also writes its bowtie2 seed index into this dir on
    first use, so a stable/writable/persistent path is exactly what we want here.

    Input:
      - db_dir:          absolute config dir (params.getorganelle_downloads)
      - force_download:  rebuild even if the sentinel exists
    Output:
      - config_dir:      the same db_dir path (val), consumed by GETORGANELLE as --config-dir
========================================================================================
*/

process DOWNLOAD_GETORGANELLE_DB {
    tag "getorganelle_db"
    label 'getorganelle_download'

    input:
    val db_dir
    val force_download

    output:
    val "${db_dir}", emit: config_dir

    script:
    """
    set -eu

    mkdir -p "${db_dir}"
    SENTINEL="${db_dir}/.getorganelle_ready"

    if [ "${force_download}" = "true" ] || [ ! -f "\${SENTINEL}" ]; then
        echo "[DOWNLOAD_GETORGANELLE_DB] Fetching all GetOrganelle databases into ${db_dir}"

        get_organelle_config.py \\
            -a all \\
            --config-dir "${db_dir}" \\
            --verbose

        # Sanity check: config dir must be non-empty after a successful fetch.
        if [ -z "\$(ls -A "${db_dir}" 2>/dev/null)" ]; then
            echo "[DOWNLOAD_GETORGANELLE_DB] ERROR: config dir empty after get_organelle_config.py" >&2
            exit 1
        fi

        date -Is > "\${SENTINEL}"
        echo "[DOWNLOAD_GETORGANELLE_DB] Done. Contents:"
        ls -la "${db_dir}"
    else
        echo "[DOWNLOAD_GETORGANELLE_DB] Databases already present at ${db_dir} (sentinel found) — skipping."
    fi
    """

    stub:
    """
    mkdir -p "${db_dir}"
    touch "${db_dir}/.getorganelle_ready"
    """
}
