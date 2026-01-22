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
  // Extract base URL from manifest URL (remove filename)
  def base_url = gxdb_manifest.replaceAll('/[^/]+$', '')
  // Extract database name (e.g., "test-only" or "all")
  def db_name = gxdb_manifest.replaceAll('.*/([^/]+)\\.manifest$', '$1')
  """
  set -euo pipefail

  mkdir -p "${gxdb_dir}"
  
  SENTINEL="${gxdb_dir}/.gxdb_ready"

  if [ "${force_download}" = "true" ] || [ ! -f "\${SENTINEL}" ] || [ -z "\$(ls -A "${gxdb_dir}" 2>/dev/null || true)" ]; then
    echo "[FCS_DB_GET] Downloading GXDB: ${db_name}"
    echo "[FCS_DB_GET] Manifest: ${gxdb_manifest}"
    echo "[FCS_DB_GET] Base URL: ${base_url}"
    echo "[FCS_DB_GET] Target directory: ${gxdb_dir}"
    
    cd "${gxdb_dir}"
    
    # Download manifest
    echo "Downloading manifest..."
    curl -L -o manifest.json "${gxdb_manifest}"
    
    # Check if jq is available for JSON parsing
    if command -v jq &> /dev/null; then
      echo "Using jq to parse manifest..."
      
      # Parse JSON manifest and extract filenames
      jq -r '.files[].fileName' manifest.json | while read -r filename; do
        if [ -n "\${filename}" ]; then
          echo "Downloading \${filename}..."
          curl -L -o "\${filename}" "${base_url}/\${filename}"
          
          if [ \$? -ne 0 ]; then
            echo "ERROR: Failed to download \${filename}"
            exit 1
          fi
          
          ls -lh "\${filename}"
        fi
      done
    else
      echo "jq not available, using fallback method..."
      
      # Fallback: Extract filenames using grep/sed for known database types
      if [ "${db_name}" = "test-only" ]; then
        FILES=(
          "test-only.gxs"
          "test-only.gxi"
          "test-only.meta.jsonl"
          "test-only.blast_div.tsv.gz"
          "test-only.seq_info.tsv.gz"
          "test-only.taxa.tsv"
        )
      elif [ "${db_name}" = "all" ]; then
        FILES=(
          "all.gxs"
          "all.gxi"
          "all.meta.jsonl"
          "all.blast_div.tsv.gz"
          "all.seq_info.tsv.gz"
          "all.taxa.tsv"
        )
      else
        echo "ERROR: Unknown database type: ${db_name}"
        echo "Please install jq or use test-only/all database"
        exit 1
      fi
      
      for filename in "\${FILES[@]}"; do
        echo "Downloading \${filename}..."
        curl -L -o "\${filename}" "${base_url}/\${filename}"
        
        if [ \$? -ne 0 ]; then
          echo "ERROR: Failed to download \${filename}"
          exit 1
        fi
        
        ls -lh "\${filename}"
      done
    fi
    
    # Verify essential files exist (.gxs and .gxi are required)
    echo "Verifying downloaded files..."
    ls -lh
    
    if ! ls *.gxs 1> /dev/null 2>&1 || ! ls *.gxi 1> /dev/null 2>&1; then
      echo "ERROR: Required database files (.gxs and .gxi) not found"
      echo "Downloaded files:"
      ls -lh
      exit 1
    fi
    
    date -Is > "\${SENTINEL}"
    echo "[FCS_DB_GET] Database download complete. Files:"
    ls -lh
  else
    echo "[FCS_DB_GET] GXDB already present at ${gxdb_dir}; skipping download."
    echo "Existing files:"
    ls -lh "${gxdb_dir}"
  fi
  """
}