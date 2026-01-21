/*
========================================================================================
    BUILD DIAMOND DATABASE MODULE
========================================================================================
    Downloads protein FASTA and taxonmap, builds DIAMOND database
    - Smart caching: skips if .dmnd already exists
    - Requires: taxdump (names.dmp, nodes.dmp)
    - Supports: gzipped or plain text inputs
========================================================================================
*/

process BUILD_DIAMOND_DB {
    tag "build_diamond_db"
    label 'diamond'
    
    input:
    path db_dir
    val db_name
    val fasta_url
    val taxonmap_url
    path taxdump_dir
    val force_build
    
    output:
    path "${db_dir}/${db_name}.dmnd", emit: dmnd
    
    script:
    """
    set -euo pipefail
    
    mkdir -p "${db_dir}"
    
    DMND="${db_dir}/${db_name}.dmnd"
    SENTINEL="${db_dir}/.${db_name}_ready"
    
    # Check if database already exists
    if [ "${force_build}" = "true" ] || \
       [ ! -s "\${DMND}" ] || \
       [ ! -f "\${SENTINEL}" ]; then
        
        echo "[BUILD_DIAMOND_DB] Building DIAMOND database: ${db_name}"
        echo "[BUILD_DIAMOND_DB] FASTA URL: ${fasta_url}"
        echo "[BUILD_DIAMOND_DB] Taxonmap URL: ${taxonmap_url}"
        
        # Download FASTA
        echo "[BUILD_DIAMOND_DB] Downloading protein sequences..."
        if [[ "${fasta_url}" == *.gz ]]; then
            curl -L "${fasta_url}" | gunzip > "${db_dir}/${db_name}.fasta"
        else
            curl -L -o "${db_dir}/${db_name}.fasta" "${fasta_url}"
        fi
        
        # Download taxonmap
        echo "[BUILD_DIAMOND_DB] Downloading taxonomy mapping..."
        if [[ "${taxonmap_url}" == *.gz ]]; then
            curl -L "${taxonmap_url}" | gunzip > "${db_dir}/${db_name}.taxonmap"
        else
            curl -L -o "${db_dir}/${db_name}.taxonmap" "${taxonmap_url}"
        fi
        
        # Build DIAMOND database
        echo "[BUILD_DIAMOND_DB] Building database (this may take several hours)..."
        diamond makedb \
            --in "${db_dir}/${db_name}.fasta" \
            --db "${db_dir}/${db_name}" \
            --taxonmap "${db_dir}/${db_name}.taxonmap" \
            --taxonnodes "${taxdump_dir}/nodes.dmp" \
            --taxonnames "${taxdump_dir}/names.dmp" \
            --threads ${task.cpus}
        
        # Cleanup intermediate files (keep only .dmnd)
        rm -f "${db_dir}/${db_name}.fasta"
        rm -f "${db_dir}/${db_name}.taxonmap"
        
        # Mark as complete
        date -Is > "\${SENTINEL}"
        
        echo "[BUILD_DIAMOND_DB] Database build complete: \${DMND}"
        ls -lh "\${DMND}"
        
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