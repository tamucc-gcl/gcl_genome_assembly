process FCS_DB_GET {
  tag "fcs_db_get"
  label 'fcs' 

  input:
    val gxdb_manifest
    val gxdb_dir
    val force_download

  output:
    val "${gxdb_dir}", emit: out_dir

  script:
  """
  set -euo pipefail

  mkdir -p "${gxdb_dir}"
  
  SENTINEL="${gxdb_dir}/.gxdb_ready"

  if [ "${force_download}" = "true" ] || [ ! -f "\${SENTINEL}" ] || [ -z "\$(ls -A "${gxdb_dir}" 2>/dev/null || true)" ]; then
    echo "[FCS_DB_GET] Downloading GXDB using singularity container"
    echo "[FCS_DB_GET] Manifest: ${gxdb_manifest}"
    echo "[FCS_DB_GET] Target directory: ${gxdb_dir}"
    
    # Nextflow automatically wraps this in singularity exec
    /app/bin/sync_files get \\
      --mft "${gxdb_manifest}" \\
      --dir "${gxdb_dir}"
    
    # Verify essential files exist
    echo "Verifying downloaded files..."
    ls -lh "${gxdb_dir}"
    
    if ! ls ${gxdb_dir}/*.gxs 1> /dev/null 2>&1 || ! ls ${gxdb_dir}/*.gxi 1> /dev/null 2>&1; then
      echo "ERROR: Required database files (.gxs and .gxi) not found"
      echo "Downloaded files:"
      ls -lh "${gxdb_dir}"
      exit 1
    fi
    
    date -Is > "\${SENTINEL}"
    echo "[FCS_DB_GET] Database download complete. Files:"
    ls -lh "${gxdb_dir}"
  else
    echo "[FCS_DB_GET] GXDB already present at ${gxdb_dir}; skipping download."
    echo "Existing files:"
    ls -lh "${gxdb_dir}"
  fi
  """
}