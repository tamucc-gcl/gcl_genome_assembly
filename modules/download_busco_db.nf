/*
========================================================================================
    DOWNLOAD BUSCO DATABASE MODULE
========================================================================================
    Downloads the specified BUSCO lineage database once, before parallel BUSCO runs.
    This prevents multiple simultaneous download attempts when running BUSCO in parallel.
    
    The downloaded database is stored in the specified output directory and can be
    reused across pipeline runs via Nextflow's caching.
========================================================================================
*/

process DOWNLOAD_BUSCO_DB {
    tag "${lineage}"
    label 'busco_download'
    
    storeDir "${params.busco_downloads}"
    
    input:
    val(lineage)
    
    output:
    path("lineages/${lineage}"), emit: db
    
    script:
    """
    # Create lineages subdirectory structure expected by BUSCO
    mkdir -p lineages
    
    # Download the specified lineage dataset
    busco --download ${lineage} --download_path .
    
    # Verify download succeeded
    if [ ! -d "lineages/${lineage}" ]; then
        echo "ERROR: Failed to download BUSCO lineage: ${lineage}"
        exit 1
    fi
    
    echo "Successfully downloaded BUSCO lineage: ${lineage}"
    """
    
    stub:
    """
    mkdir -p lineages/${lineage}
    touch lineages/${lineage}/dataset.cfg
    """
}