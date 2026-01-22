/*
========================================================================================
    BUILD DIAMOND DATABASE MODULE
========================================================================================
    Downloads protein FASTA and taxonmap, builds DIAMOND database
    - Smart caching: skips if .dmnd already exists
    - Requires: taxdump (names.dmp, nodes.dmp)
    - Supports: gzipped or plain text inputs
    - Works like FCS_DB_GET: manages database in permanent location
========================================================================================
*/

process BUILD_DIAMOND_DB {
    tag "build_diamond_db"
    label 'diamond'
    
    input:
    val db_dir        // Directory path as a string
    val db_name
    val fasta_url
    val taxonmap_url
    path taxdump_dir  // This gets staged, but that's OK since it's small
    val force_build
    
    output:
    val "${db_dir}/${db_name}.dmnd", emit: dmnd  // Output path as value, not file
    
    script:
    """
    set -euo pipefail
    set -o pipefail  # Ensure pipe failures are caught
    
    # Ensure target directory exists
    mkdir -p "${db_dir}"
    
    DMND="${db_dir}/${db_name}.dmnd"
    SENTINEL="${db_dir}/.${db_name}_ready"
    
    # Check if database already exists
    if [ "${force_build}" = "true" ] || \
       [ ! -s "\${DMND}" ] || \
       [ ! -f "\${SENTINEL}" ]; then
        
        echo "[BUILD_DIAMOND_DB] Building DIAMOND database: ${db_name}"
        echo "[BUILD_DIAMOND_DB] Target location: ${db_dir}"
        echo "[BUILD_DIAMOND_DB] FASTA URL: ${fasta_url}"
        echo "[BUILD_DIAMOND_DB] Taxonmap URL: ${taxonmap_url}"
        echo ""
        echo "WARNING: This will download ~80GB (compressed) and uncompress to ~250GB"
        echo "         Final database will be ~120GB. Ensure sufficient disk space!"
        echo ""
        
        # Download FASTA with better error handling
        echo "[BUILD_DIAMOND_DB] Downloading protein sequences (this will take hours)..."
        echo "[BUILD_DIAMOND_DB] Started: \$(date)"
        
        if [[ "${fasta_url}" == *.gz ]]; then
            # Use -f to fail on HTTP errors, --progress-bar for visibility
            if ! curl -f -L --progress-bar "${fasta_url}" | gunzip > "${db_dir}/${db_name}.fasta"; then
                echo "[BUILD_DIAMOND_DB] ERROR: FASTA download failed!" >&2
                rm -f "${db_dir}/${db_name}.fasta"
                exit 1
            fi
        else
            if ! curl -f -L --progress-bar -o "${db_dir}/${db_name}.fasta" "${fasta_url}"; then
                echo "[BUILD_DIAMOND_DB] ERROR: FASTA download failed!" >&2
                rm -f "${db_dir}/${db_name}.fasta"
                exit 1
            fi
        fi
        
        # Check downloaded file size
        FASTA_SIZE=\$(stat -c%s "${db_dir}/${db_name}.fasta" 2>/dev/null || stat -f%z "${db_dir}/${db_name}.fasta" 2>/dev/null || echo 0)
        FASTA_SIZE_GB=\$(echo "scale=2; \${FASTA_SIZE}/1024/1024/1024" | bc)
        echo "[BUILD_DIAMOND_DB] Downloaded FASTA size: \${FASTA_SIZE_GB} GB"
        
        # Sanity check - nr.fasta should be at least 200 GB
        if [ "\${FASTA_SIZE}" -lt 200000000000 ]; then
            echo "[BUILD_DIAMOND_DB] ERROR: FASTA file too small (\${FASTA_SIZE_GB} GB)!" >&2
            echo "[BUILD_DIAMOND_DB] Expected at least 200 GB for NCBI nr" >&2
            echo "[BUILD_DIAMOND_DB] Download may have been interrupted" >&2
            rm -f "${db_dir}/${db_name}.fasta"
            exit 1
        fi
        
        echo "[BUILD_DIAMOND_DB] FASTA download complete: \$(date)"
        
        # Download taxonmap with better error handling
        echo "[BUILD_DIAMOND_DB] Downloading taxonomy mapping..."
        echo "[BUILD_DIAMOND_DB] Started: \$(date)"
        
        if [[ "${taxonmap_url}" == *.gz ]]; then
            if ! curl -f -L --progress-bar "${taxonmap_url}" | gunzip > "${db_dir}/${db_name}.taxonmap"; then
                echo "[BUILD_DIAMOND_DB] ERROR: Taxonmap download failed!" >&2
                rm -f "${db_dir}/${db_name}.taxonmap"
                exit 1
            fi
        else
            if ! curl -f -L --progress-bar -o "${db_dir}/${db_name}.taxonmap" "${taxonmap_url}"; then
                echo "[BUILD_DIAMOND_DB] ERROR: Taxonmap download failed!" >&2
                rm -f "${db_dir}/${db_name}.taxonmap"
                exit 1
            fi
        fi
        
        TAXONMAP_SIZE=\$(stat -c%s "${db_dir}/${db_name}.taxonmap" 2>/dev/null || stat -f%z "${db_dir}/${db_name}.taxonmap" 2>/dev/null || echo 0)
        TAXONMAP_SIZE_GB=\$(echo "scale=2; \${TAXONMAP_SIZE}/1024/1024/1024" | bc)
        echo "[BUILD_DIAMOND_DB] Downloaded taxonmap size: \${TAXONMAP_SIZE_GB} GB"
        echo "[BUILD_DIAMOND_DB] Taxonmap download complete: \$(date)"
        
        # Build DIAMOND database directly in target directory
        echo ""
        echo "[BUILD_DIAMOND_DB] Building DIAMOND database (this will take many hours)..."
        echo "[BUILD_DIAMOND_DB] Started: \$(date)"
        echo "[BUILD_DIAMOND_DB] Using ${task.cpus} threads"
        echo ""
        
        if ! diamond makedb \
            --in "${db_dir}/${db_name}.fasta" \
            --db "${db_dir}/${db_name}" \
            --taxonmap "${db_dir}/${db_name}.taxonmap" \
            --taxonnodes "${taxdump_dir}/nodes.dmp" \
            --taxonnames "${taxdump_dir}/names.dmp" \
            --threads ${task.cpus}; then
            echo "[BUILD_DIAMOND_DB] ERROR: DIAMOND makedb failed!" >&2
            exit 1
        fi
        
        echo "[BUILD_DIAMOND_DB] DIAMOND makedb complete: \$(date)"
        
        # Cleanup intermediate files (keep only .dmnd)
        echo "[BUILD_DIAMOND_DB] Cleaning up intermediate files..."
        rm -f "${db_dir}/${db_name}.fasta"
        rm -f "${db_dir}/${db_name}.taxonmap"
        
        # Mark as complete
        date -Is > "\${SENTINEL}"
        
        echo ""
        echo "[BUILD_DIAMOND_DB] Database build complete: \${DMND}"
        ls -lh "\${DMND}"
        echo ""
        
    else
        echo "[BUILD_DIAMOND_DB] DIAMOND database already exists: \${DMND}"
        echo "[BUILD_DIAMOND_DB] Skipping build"
        ls -lh "\${DMND}"
    fi
    """
    
    stub:
    """
    mkdir -p "${db_dir}"
    touch "${db_dir}/${db_name}.dmnd"
    touch "${db_dir}/.${db_name}_ready"
    """
}