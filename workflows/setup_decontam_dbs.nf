/*
========================================================================================
    SETUP DECONTAMINATION DATABASES WORKFLOW
========================================================================================
    Purpose:
    - Downloads/prepares FCS-GX database once for all samples
    - Downloads/prepares DIAMOND database if evidence is enabled
    - Downloads NCBI taxonomy database
    - Runs in parallel with BAM conversion and assembly
    
    Design:
    - Single execution at pipeline start
    - Shared databases for all downstream decontamination
    - Smart caching: skip downloads if databases already exist
    - Clear separation: taxdump is independent of DIAMOND database
========================================================================================
*/

nextflow.enable.dsl=2

include { FCS_DB_GET }         from '../modules/fcs_dbgx.nf'
include { DOWNLOAD_TAXDUMP }   from '../modules/download_taxdump.nf'
include { BUILD_DIAMOND_DB }   from '../modules/build_diamond_db.nf'

workflow SETUP_DECONTAM_DBS {
    main:
    /*
    // DEBUG
    log.info """
    ================================================================================
    SETUP_DECONTAM_DBS DEBUG
    ================================================================================
    params.decon.make_blobtools_evidence: ${params.decon?.make_blobtools_evidence}
    Check result: ${params.decon?.make_blobtools_evidence ?: false}
    params.diamond.dmnd: ${params.diamond?.dmnd}
    ================================================================================
    """
    */

    /*
    ========================================================================================
        FCS-GX Database Setup
    ========================================================================================
    */
    def gxdb_profile  = params.gxdb?.profile  ?: 'all'
    def gxdb_manifest = params.gxdb?.manifest ?: (
        gxdb_profile == 'test-only'
            ? 'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest'
            : 'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest'
    )
    def gxdb_dir_path = file(params.gxdb?.dir ?: "${params.db_base}/fcs-gx")
    def gxdb_force = (params.gxdb?.force ?: false) as boolean

    // Create channels for the inputs
    ch_gxdb_manifest = Channel.value(gxdb_manifest)
    ch_gxdb_dir = Channel.value(gxdb_dir_path.toString())  // ← Convert to string
    ch_gxdb_force = Channel.value(gxdb_force)

    // Call the process
    FCS_DB_GET(
        ch_gxdb_manifest,
        ch_gxdb_dir,
        ch_gxdb_force
    )

    /*
    ========================================================================================
        Optional: Evidence Databases (DIAMOND + Taxonomy)
    ========================================================================================
    */
    if (params.decon?.make_blobtools_evidence ?: false) {
        
        def taxdump_dir_path = file(params.diamond?.taxdump_dir ?: './db/taxdump')
        def taxdump_force = (params.diamond?.force ?: false) as boolean
        
        // Create channels
        ch_taxdump_dir = Channel.value(taxdump_dir_path)
        ch_taxdump_force = Channel.value(taxdump_force)
        
        // Always ensure taxdump exists (needed by both DIAMOND and BlobTools)
        DOWNLOAD_TAXDUMP(
            ch_taxdump_dir,
            ch_taxdump_force
        )
        
        // Check if user provided pre-built DIAMOND database
        def prebuilt_dmnd = params.diamond?.dmnd ? file(params.diamond.dmnd) : null
        
        if (prebuilt_dmnd && prebuilt_dmnd.exists()) {
            // Use pre-built database - no download/build needed
            diamond_db_out = Channel.value(prebuilt_dmnd)
            
        } else {
            // Need to build DIAMOND database
            def diamond_dir = file(params.diamond?.dir ?: './db/diamond')
            def diamond_name = params.diamond?.name ?: 'proteins'
            def fasta_url = params.diamond?.fasta_url
            def taxonmap_url = params.diamond?.taxonmap_url
            def diamond_force = (params.diamond?.force ?: false) as boolean
            
            // Validate that URLs are provided
            if (!fasta_url || !taxonmap_url) {
                error """
                ================================================================================
                ERROR: DIAMOND database build requested but URLs not provided!
                
                You must either:
                1. Provide a pre-built database:
                   --diamond.dmnd /path/to/existing.dmnd
                
                2. Provide URLs to build a new database:
                   --diamond.fasta_url https://url/to/proteins.fasta.gz
                   --diamond.taxonmap_url https://url/to/taxonmap.gz
                
                Example (NCBI nr):
                   --diamond.fasta_url https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
                   --diamond.taxonmap_url https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
                ================================================================================
                """
            }
            
            BUILD_DIAMOND_DB(
                Channel.value(diamond_dir),
                Channel.value(diamond_name),
                Channel.value(fasta_url),
                Channel.value(taxonmap_url),
                DOWNLOAD_TAXDUMP.out.taxdump_dir,
                Channel.value(diamond_force)
            )
            
            diamond_db_out = BUILD_DIAMOND_DB.out.dmnd
        }
        
        taxdump_out = DOWNLOAD_TAXDUMP.out.taxdump_dir
        
    } else {
        // Evidence not requested - emit empty channels
        diamond_db_out = Channel.empty()
        taxdump_out = Channel.empty()
    }

    emit:
    gxdb_dir = FCS_DB_GET.out.out_dir
    diamond_db = diamond_db_out
    taxdump_dir = taxdump_out
}