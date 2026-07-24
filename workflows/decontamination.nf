nextflow.enable.dsl=2

/*
 * Core FCS decontamination
 */
include { FCS_DB_GET }          from '../modules/fcs_dbgx.nf'
include { FCS_ADAPTOR }         from '../modules/fcs_adaptor.nf'
include { FCS_GX_SCREEN }       from '../modules/fcs_gx_screen.nf'
include { FCS_CLEAN_GENOME }    from '../modules/fcs_clean_genome.nf'

/*
 * Optional evidence (coverage + taxonomy + blobtools plots)
 */
include { MAP_READS_MINIMAP2 }  from '../modules/map_reads_minimap2.nf'
include { DIAMOND_DB_GET }      from '../modules/diamond_db_get.nf'
include { DIAMOND_BLASTX }      from '../modules/diamond_blastx.nf'
include { BLOBTOOLS_CREATE }    from '../modules/blobtools2_create.nf'
include { BLOBTOOLS_VIEWPLOT }  from '../modules/blobtools2_viewplot.nf'

workflow DECONTAM_FCS_AUTO {

  take:
    assembly_fa
    hifi_reads

  main:

    /*
     * Resolve GXDB manifest (can be overridden explicitly in params.gxdb_manifest)
     */
    def gxdb_profile  = params.gxdb_profile  ?: 'all'        // 'all' | 'test-only'
    def gxdb_manifest = params.gxdb_manifest ?: (
      gxdb_profile == 'test-only'
        ? 'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest'
        : 'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest'
    )
    def gxdb_dir   = file(params.gxdb_dir ?: './db/fcs-gx')
    def gxdb_force = (params.gxdb_force ?: false) as boolean

    /*
     * Ensure GXDB exists (download if missing/empty or force=true)
     */
    gxdb = FCS_DB_GET(
      gxdb_manifest,
      gxdb_dir,
      gxdb_force
    )

    /*
     * Optional adaptor/vector screening
     */
    cleaned_input = assembly_fa
    if (params.run_fcs_adaptor ?: false) {
      ad = FCS_ADAPTOR(
        assembly_fa,
        (params.decon_fcsadaptor_mode ?: 'euk'),
        (params.decon_container_engine ?: 'singularity')
      )
      cleaned_input = ad.cleaned_fasta
    }

    /*
     * FCS-GX screen + clean
     */
    gx = FCS_GX_SCREEN(
      cleaned_input,
      (params.decon_source_taxid ?: 7898),  // 7898 = Actinopterygii fallback
      gxdb.out_dir,
      (params.decon?.cpus ?: 32)
    )

    clean = FCS_CLEAN_GENOME(
      cleaned_input,
      gx.action_report
    )

    /*
     * Optional evidence pack (blobtools2) — runs on cleaned assembly
     */
    if (params.run_blobtools_evidence ?: false) {

      // ---- DIAMOND DB resolution ----
      // Preferred: user provides prebuilt .dmnd at params.diamond_dmnd
      // Else: download/build to params.diamond_dir (and ensure taxdump exists)
      diamond_db_path = params.diamond_dmnd ? file(params.diamond_dmnd) : null
      taxdump_dir_path = file(params.diamond_taxdump_dir ?: './db/taxdump')

      if ( !diamond_db_path ) {
        def diamond_profile = params.diamond_profile ?: 'custom'  // 'custom' is safest unless you pin URLs
        def diamond_dir     = file(params.diamond_dir ?: './db/diamond')
        def diamond_name    = params.diamond_name ?: 'proteins'

        diamond_prep = DIAMOND_DB_GET(
          diamond_profile,
          diamond_dir,
          diamond_name,
          params.diamond_fasta_url,
          params.diamond_taxonmap_url,
          taxdump_dir_path,
          (params.diamond_force ?: false) as boolean,
          (params.decon?.cpus ?: 32)
        )

        diamond_db_path  = diamond_prep.dmnd
        taxdump_dir_path = diamond_prep.taxdump_out
      }

      // ---- Coverage from HiFi (evidence only) ----
      map = MAP_READS_MINIMAP2(
        clean.decontaminated_fasta,
        hifi_reads,
        (params.evidence_map_preset ?: 'map-hifi'),
        (params.decon?.cpus ?: 32)
      )

      // ---- DIAMOND blastx (evidence taxonomy signal) ----
      hits = DIAMOND_BLASTX(
        clean.decontaminated_fasta,
        diamond_db_path,
        (params.decon?.cpus ?: 32),
        (params.evidence_diamond_max_target_seqs ?: 1),
        (params.evidence_diamond_evalue ?: 1e-25)
      )

      // ---- blobtools2 evidence outputs ----
      blob = BLOBTOOLS_CREATE(
        clean.decontaminated_fasta,
        hits.out_hits,
        map.bam,
        taxdump_dir_path,
        (params.evidence_blob_min_contig_bp ?: 1000)
      )

      viz = BLOBTOOLS_VIEWPLOT(blob.out_blobdir)
    }

  emit:
    decontaminated_fasta = clean.decontaminated_fasta
    contaminants_fasta   = clean.contaminants_fasta
    action_report        = gx.action_report
    taxonomy_report      = gx.taxonomy_report
    fcs_stdout_summary   = gx.stdout_log

    // Optional evidence outputs
    coverage_bam         = map.bam          optional true
    blobtools_blobdir    = blob.out_blobdir optional true
    blobtools_outputs    = viz.out_dir      optional true
}
