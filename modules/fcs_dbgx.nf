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

  # Resolve the actual path (follow symlinks)
  ACTUAL_DIR=\$(readlink -f "${gxdb_dir}" 2>/dev/null || echo "${gxdb_dir}")
  
  # Create the actual directory
  mkdir -p "\${ACTUAL_DIR}"
  
  SENTINEL="\${ACTUAL_DIR}/.gxdb_ready"

  if [ "${force_download}" = "true" ] || [ ! -f "\${SENTINEL}" ] || [ -z "\$(ls -A "\${ACTUAL_DIR}" 2>/dev/null || true)" ]; then
    echo "[FCS_DB_GET] Downloading GXDB using manifest: ${gxdb_manifest}"
    echo "[FCS_DB_GET] Target directory: \${ACTUAL_DIR}"
    
    cd "\${ACTUAL_DIR}"
    
    # Download manifest file
    echo "Downloading manifest..."
    curl -L -o manifest.txt "${gxdb_manifest}"
    
    # Parse manifest and download each file
    echo "Parsing manifest and downloading files..."
    while IFS=\$'\\t' read -r filename url hash; do
      # Skip comments and empty lines
      if [[ "\${filename}" =~ ^#.*\$ ]] || [[ -z "\${filename}" ]]; then
        continue
      fi
      
      if [ -n "\${url}" ]; then
        echo "Downloading \${filename} from \${url}..."
        curl -L -o "\${filename}" "\${url}"
        
        if [ \$? -ne 0 ]; then
          echo "ERROR: Failed to download \${filename}"
          exit 1
        fi
      fi
    done < manifest.txt
    
    # Verify essential files exist
    echo "Verifying downloaded files..."
    ls -lh
    
    # Check for expected database files
    file_count=\$(ls *.gx* 2>/dev/null | wc -l)
    if [ "\${file_count}" -eq 0 ]; then
      echo "ERROR: No .gx* database files found after download"
      exit 1
    fi
    
    date -Is > "\${SENTINEL}"
    echo "[FCS_DB_GET] Database download complete. Files:"
    ls -lh *.gx* 2>/dev/null || ls -lh
  else
    echo "[FCS_DB_GET] GXDB already present at \${ACTUAL_DIR}; skipping download."
    echo "Existing files:"
    ls -lh "\${ACTUAL_DIR}" 2>/dev/null || echo "Directory is empty or inaccessible"
  fi
  """
}