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
  // Extract database name from manifest (test-only or all)
  def db_name = gxdb_manifest.replaceAll('.*/([^/]+)\\.manifest.*$', '$1')
  """
  set -euo pipefail

  mkdir -p "${gxdb_dir}"
  
  # Make sentinel specific to database type
  SENTINEL="${gxdb_dir}/.gxdb_ready_${db_name}"

  if [ "${force_download}" = "true" ] || [ ! -f "\${SENTINEL}" ] || [ ! ls ${gxdb_dir}/${db_name}.gx* 1> /dev/null 2>&1 ]; then
    echo "[FCS_DB_GET] Downloading GXDB: ${db_name}"
    echo "[FCS_DB_GET] Manifest: ${gxdb_manifest}"
    echo "[FCS_DB_GET] Target directory: ${gxdb_dir}"
    
    # Nextflow automatically wraps this in singularity exec
    /app/bin/sync_files get \\
      --mft "${gxdb_manifest}" \\
      --dir "${gxdb_dir}"
    
    # Verify essential files exist for this specific database
    echo "Verifying downloaded files..."
    ls -lh "${gxdb_dir}"
    
    if ! ls ${gxdb_dir}/${db_name}.gxs 1> /dev/null 2>&1 || ! ls ${gxdb_dir}/${db_name}.gxi 1> /dev/null 2>&1; then
      echo "ERROR: Required database files (${db_name}.gxs and ${db_name}.gxi) not found"
      echo "Downloaded files:"
      ls -lh "${gxdb_dir}"
      exit 1
    fi
    
    date -Is > "\${SENTINEL}"
    echo "[FCS_DB_GET] Database download complete. Files:"
    ls -lh "${gxdb_dir}"/${db_name}*
  else
    echo "[FCS_DB_GET] GXDB '${db_name}' already present at ${gxdb_dir}; skipping download."
    echo "Existing files:"
    ls -lh "${gxdb_dir}"/${db_name}* 2>/dev/null || echo "No files found for ${db_name}"
  fi
  """
}