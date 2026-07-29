/*
========================================================================================
    DOWNLOAD MITOS2 REFERENCE DATA MODULE
========================================================================================
    Repo location: modules/download_mitos_db.nf

    Fetches the MITOS2 RefSeq reference set(s) once into a permanent shared dir, before any
    parallel ANNOTATE_MITO runs. Mirrors DOWNLOAD_BUSCO_DB: download into an absolute dir,
    guard with a sentinel, skip if present, emit the parent dir for `runmitos.py -R`.

    Source (verified): Zenodo record 3685310 ("Reference data for MITOS2", incl. refseq89).
      refseq89m.tar.bz2  Metazoa (animal mito)   ~20 MB
      refseq89f.tar.bz2  Fungi                    ~2 MB
      refseq89o.tar.bz2  Opisthokonta            ~21 MB   (also refseq63/refseq39 available)
    Each tarball extracts to a directory named after the set (e.g. refseq89m/), so the
    download dir becomes the `-R` parent and the set name is the `-r` argument.

    Input:
      - db_dir:        absolute parent dir (params.mitos_downloads)
      - refseq_sets:   comma-separated set names to fetch (params.mitos_refseq_sets)
      - force:         re-download even if the sentinel exists
    Output:
      - db:            db_dir (val), consumed by ANNOTATE_MITO as -R
========================================================================================
*/

process DOWNLOAD_MITOS_DB {
    tag "mitos_db"
    label 'mitos_download'

    input:
    val db_dir
    val refseq_sets
    val force

    output:
    val "${db_dir}", emit: db

    script:
    def base = 'https://zenodo.org/records/3685310/files'
    """
    set -eu
    mkdir -p "${db_dir}"
    SENTINEL="${db_dir}/.mitos_ready"

    if [ "${force}" = "true" ] || [ ! -f "\${SENTINEL}" ]; then
        echo "${refseq_sets}" | tr ',' '\\n' | while read -r s; do
            s=\$(echo "\$s" | tr -d '[:space:]')
            [ -z "\$s" ] && continue
            if [ ! -d "${db_dir}/\$s" ]; then
                echo "[DOWNLOAD_MITOS_DB] fetching \$s from Zenodo"
                curl -fsSL -o "${db_dir}/\$s.tar.bz2" "${base}/\$s.tar.bz2?download=1"
                tar -xjf "${db_dir}/\$s.tar.bz2" -C "${db_dir}"
                rm -f "${db_dir}/\$s.tar.bz2"
            fi
            [ -d "${db_dir}/\$s" ] || echo "[DOWNLOAD_MITOS_DB] WARNING: ${db_dir}/\$s missing after extract"
        done
        date -Is > "\${SENTINEL}"
        echo "[DOWNLOAD_MITOS_DB] done:"; ls -la "${db_dir}"
    else
        echo "[DOWNLOAD_MITOS_DB] present (sentinel found) — skipping"
    fi
    """

    stub:
    """
    mkdir -p "${db_dir}/refseq89m"
    touch "${db_dir}/.mitos_ready"
    """
}
