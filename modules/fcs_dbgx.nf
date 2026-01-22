process FCS_DB_GET {
  tag "fcs_db_get"
  label 'fcs' 

  input:
    val gxdb_manifest
    path gxdb_dir
    val force_download

  output:
    path "${gxdb_dir}", emit: out_dir

  script:
  """
  set -euo pipefail

  command -v run_gx.py >/dev/null 2>&1 || { echo "ERROR: run_gx.py not found in PATH"; exit 127; }

  mkdir -p "${gxdb_dir}"

  # Consider DB "present" if directory has files AND sentinel exists
  SENTINEL="${gxdb_dir}/.gxdb_ready"

  if [ "${force_download}" = "true" ] || [ ! -f "\${SENTINEL}" ] || [ -z "\$(ls -A "${gxdb_dir}" 2>/dev/null || true)" ]; then
    echo "[FCS_DB_GET] Downloading GXDB using manifest: ${gxdb_manifest}"
    run_gx.py db get --mft "${gxdb_manifest}" --dir "${gxdb_dir}"
    date -Is > "\${SENTINEL}"
  else
    echo "[FCS_DB_GET] GXDB already present at ${gxdb_dir}; skipping download."
  fi
  """
}
