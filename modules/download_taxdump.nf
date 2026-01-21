/*
========================================================================================
    DOWNLOAD TAXDUMP MODULE
========================================================================================
    Downloads NCBI taxonomy database (taxdump)
    - Required by: DIAMOND makedb, BlobTools2
    - Contains: names.dmp, nodes.dmp (taxonomic names and tree structure)
    - Smart caching: skips download if files already exist
========================================================================================
*/

process DOWNLOAD_TAXDUMP {
    tag "download_taxdump"
    label 'download_db'
    
    input:
    path taxdump_dir
    val force_download
    
    output:
    path "${taxdump_dir}", emit: taxdump_dir
    
    script:
    """
    set -euo pipefail
    
    mkdir -p "${taxdump_dir}"
    
    SENTINEL="${taxdump_dir}/.taxdump_ready"
    
    # Check if taxdump already exists and is complete
    if [ "${force_download}" = "true" ] || \
       [ ! -s "${taxdump_dir}/names.dmp" ] || \
       [ ! -s "${taxdump_dir}/nodes.dmp" ] || \
       [ ! -f "\${SENTINEL}" ]; then
        
        echo "[DOWNLOAD_TAXDUMP] Downloading NCBI taxdump to ${taxdump_dir}"
        
        # Download taxdump from NCBI FTP
        curl -L -o "${taxdump_dir}/taxdump.tar.gz" \
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        
        # Extract only the files we need
        tar -xzf "${taxdump_dir}/taxdump.tar.gz" \
            -C "${taxdump_dir}" \
            names.dmp nodes.dmp
        
        # Cleanup
        rm -f "${taxdump_dir}/taxdump.tar.gz"
        
        # Mark as complete
        date -Is > "\${SENTINEL}"
        
        echo "[DOWNLOAD_TAXDUMP] Download complete"
    else
        echo "[DOWNLOAD_TAXDUMP] Taxdump already present at ${taxdump_dir} - skipping download"
        echo "[DOWNLOAD_TAXDUMP] Found:"
        echo "  - names.dmp: \$(wc -l < "${taxdump_dir}/names.dmp") lines"
        echo "  - nodes.dmp: \$(wc -l < "${taxdump_dir}/nodes.dmp") lines"
    fi
    """
    
    stub:
    """
    mkdir -p "${taxdump_dir}"
    touch "${taxdump_dir}/names.dmp"
    touch "${taxdump_dir}/nodes.dmp"
    touch "${taxdump_dir}/.taxdump_ready"
    """
}