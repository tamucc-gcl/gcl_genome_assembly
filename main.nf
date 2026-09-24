#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    Genome Assembly and Scaffolding Pipeline
========================================================================================
    Author: Jason Selwyn
    Description: Pipeline for assembling and scaffolding genomes from HiFi and Hi-C reads
    
    MODULAR DESIGN:
    - Maximum parallelization
    - Reusable QC workflows as "functions"
    - Clear separation of concerns
    - Optional misassembly correction (Inspector)
    - Optional decontamination at multiple stages
    - Conditional iterative scaffolding (second round only when beneficial)
    - Gap filling using TGSGapCloser
========================================================================================
*/

/*
========================================================================================
    PARAMETERS - FLATTENED FOR EASY COMMAND-LINE OVERRIDE
========================================================================================
    All parameters can now be overridden individually without affecting others.
    Example: --decon_run_on_contigs true --decon_source_taxid 373251
    
    IMPORTANT: Place these BEFORE any conditional logic that uses params
========================================================================================
*/

// ============================================================================
// Core pipeline parameters (keep existing ones)
// ============================================================================

// Shared databases are snapshots. storeDir reuses existing outputs even when the
// task script or source URL changes; select a new directory to refresh a database.
if (params.diamond_force || params.gxdb_force || params.getorganelle_force_download || params.mitos_force_download) {
    error 'Database force-download flags are retired. Set the relevant database directory to a new snapshot directory instead.'
}
// Print pipeline header
// helpers for the startup banner
def onOff  = { it ? 'Enabled' : 'Disabled' }
def stages = { c, s -> ([c ? 'Contigs' : null, s ? 'Scaffolds' : null].findAll().join(' + ')) ?: 'Disabled' }
def orAuto = { v, src -> (v == null || v == 'auto') ? "auto (${src})" : v.toString() }

log.info """\
    ============================================================
     GENOME ASSEMBLY PIPELINE  ${workflow.manifest.version ?: ''}
    ============================================================
     Run              : ${workflow.runName}${workflow.revision ? "  (${workflow.revision})" : ''}
     Profile          : ${workflow.profile}
     Sample sheet     : ${params.sample_sheet}
     Output dir       : ${params.outdir}
    ------------------------------------------------------------
     purge_dups       : ${onOff(params.run_purge_dups)}
     Pilon polish     : ${onOff(params.run_pilon)}
     Inspector        : ${stages(params.run_inspector_contigs, params.run_inspector_scaffolds)}
     Decontamination  : ${stages(params.run_decon_contigs, params.run_decon_scaffolds)}${params.run_fcs_adaptor ? '  (+FCS-adaptor)' : ''}
     Scaffold round 2 : ${onOff(params.run_scaffold_round2)}
     Gap filling      : Enabled
     Hi-C contact maps: ${onOff(params.run_post_assembly != false && params.run_final_contact_maps)}
     Pairwise synteny : ${params.run_post_assembly != false && params.run_pairwise_alignments ? params.pairwise_alignment_mode : 'Disabled'}
     Teloclip extend  : ${onOff(params.run_teloclip_extend)}
    ------------------------------------------------------------
     Analysis k-mer   : ${params.kmer_size}
     Genetic code     : ${orAuto(params.mitochondria_genetic_code , 'from taxid')}
     Telomere motif   : ${orAuto(params.telomere_motif, 'from taxid')}
     Haploid genome size: ${orAuto(params.haploid_genome_size, 'from estimate')}
    ============================================================
    """.stripIndent()

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!params.sample_sheet) {
    exit 1, "Sample sheet not specified! Please provide --sample_sheet"
}
if (params.run_blobtools_evidence) {
    error 'BlobTools/decontamination evidence is unsupported pending its rebuild; set --run_blobtools_evidence false'
}

/*
========================================================================================
    IMPORT FUNCTIONS
========================================================================================
*/
include { parseSampleSheet } from './functions/input_validation.nf'
include { INPUT_PREPARATION } from './workflows/input_preparation.nf'
include { ASSEMBLY_RUN_SUMMARY } from './modules/assembly_run_summary.nf'
include { forkHaplotypeMeta } from './functions/meta.nf'

/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/
include { HIC_QC as HIC_QC_RAW } from './workflows/hic_qc.nf'
include { HIC_QC as HIC_QC_TRIMMED } from './workflows/hic_qc.nf'
include { HIFI_QC } from './workflows/hifi_qc.nf'

// Assembly
include { HIC_SCAFFOLDING } from './workflows/hic_scaffolding.nf'
include { READ_PREPARATION } from './workflows/read_preparation.nf'
include { CONTIG_REFINEMENT } from './workflows/contig_refinement.nf'
include { CONTIG_ASSEMBLY } from './workflows/contig_assembly.nf'

//Organelle assembly/annotation
include { ORGANELLE } from './workflows/organelle.nf'

//Short-read qc #Hi-C, SSL, TellSeq
include { SHORTREAD_QC as SHORTREAD_QC_RAW }     from './workflows/shortread_qc.nf'
include { SHORTREAD_QC as SHORTREAD_QC_TRIMMED } from './workflows/shortread_qc.nf'

// Assembly QC
include { QC_PHASE } from './workflows/qc_phase.nf'

// HI-C MODULAR WORKFLOWS
include { HIC_QC_FROM_BAM as HIC_QC_FROM_BAM_RAW } from './workflows/hic_qc_from_bam.nf'
include { HIC_QC_FROM_BAM as HIC_QC_FROM_BAM_FILTERED } from './workflows/hic_qc_from_bam.nf'
include { HIC_SCAFFOLD_QC } from './workflows/hic_scaffold_qc.nf'

// DECONTAMINATION MODULAR WORKFLOWS
include { SETUP_DECONTAM_DBS } from './workflows/setup_decontam_dbs.nf'


//Final Visualization
include { FINAL_VIZ } from './workflows/final_viz.nf'
include { FINAL_HIC_MAPS } from './workflows/final_hic_maps.nf'

// Pangenome
include { PANGENOME } from './workflows/pangenome.nf'

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { RESOLVE_TAXONOMY } from './modules/resolve_taxonomy.nf'
include { buscoLineageFor; kingdomFlag; organismName; geneticCodeFor; telomereMotifFor; organelleTypesFor; getorganelleRecursionFor; getorganelleKmersFor; getorganelleCoverageFor } from './functions/taxonomy.nf'
include { DOWNLOAD_TAXDUMP } from './modules/download_taxdump.nf'

include { BUILD_MERYL_DB } from './modules/build_meryl_db.nf'

include { FIND_MITO_REFERENCE } from './modules/find_mito_reference.nf'


include { SNAIL_PLOT as SNAIL_PLOT_FINAL } from './modules/snail_plot.nf'
//include { SCAN_TELOMERES; COLLECT_TELOMERE_RESULTS } from './modules/scan_telomeres.nf'
include { DOWNLOAD_BUSCO_DB } from './modules/download_busco_db.nf'
include { COVERAGE_BOOK } from './modules/coverage_book.nf'
include { REPORTING } from './workflows/reporting.nf'
include { FINALIZE_ASSEMBLY } from './modules/finalize_assembly.nf'
include { HARMONIZE_SCAFFOLDS } from './workflows/harmonize_scaffolds.nf'
include { CHIMERA } from './workflows/chimera.nf'
include { COLLECT_NAME_MAPS } from './modules/collect_name_maps.nf'

// ── helper scripts declared as inputs so edits invalidate the cache ──
ch_compile_qc_script      = file("${projectDir}/r_scripts/compile_qc.R",                         checkIfExists: true)
ch_summary_report_script  = file("${projectDir}/r_scripts/generate_summary_report.R",            checkIfExists: true)
ch_dotplot_script         = file("${projectDir}/r_scripts/dotplot_paf.R",                        checkIfExists: true)
ch_riparian_script        = file("${projectDir}/r_scripts/riparian_paf.R",                       checkIfExists: true)
ch_assembly_report_script = file("${projectDir}/py_scripts/generate_assembly_report.py",         checkIfExists: true)
ch_coverage_book_script   = file("${projectDir}/py_scripts/bigwig_genome_book.py",               checkIfExists: true)
ch_tad_book_script        = file("${projectDir}/py_scripts/make_tad_book.py",                    checkIfExists: true)
ch_compartments_script    = file("${projectDir}/py_scripts/plot_compartments_pc1_genomewide.py", checkIfExists: true)
ch_harmonize_script       = file("${projectDir}/py_scripts/harmonize_names.py",                     checkIfExists: true)
/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

// User-run input checkpoint: performs input validation/reporting only.
workflow VALIDATE_INPUTS {
    def validated_inputs = parseSampleSheet(params.sample_sheet, params.hic_readsets)
    def capabilities = validated_inputs.capabilities
    log.info '[INPUT] Enabled branches: ' + capabilities.findAll { key, enabled -> enabled }.keySet().join(', ')
    INPUT_PREPARATION(validated_inputs)
}

workflow {
    
    // ---- run-scope gates, defined once before any use ----
    def qc_mode     = params.qc_mode ?: 'all_stages'
    def qc_on       = ( qc_mode != 'none' )
    def run_all_qc  = qc_on && ( qc_mode != 'final_only' )
    def post_on     = ( params.run_post_assembly != false )
    def chimera_on  = params.chimera_break && params.chimera_break.toString() != 'false'

    // Parse sample sheet -> per-sample tuple(meta, reads)
    def validated_inputs = parseSampleSheet(params.sample_sheet, params.hic_readsets)
    def capabilities = validated_inputs.capabilities
    log.info '[INPUT] Enabled branches: ' + capabilities.findAll { key, enabled -> enabled }.keySet().join(', ')
    INPUT_PREPARATION(validated_inputs)
    ch_input = INPUT_PREPARATION.out.samples
    ch_input_status = INPUT_PREPARATION.out.status

    // ploidy + haploid-size overrides ride sample-keyed side-channels (genomescope/hifiasm only),
    // so ploidy/size tweaks stay off the per-sample task hash.
    ch_ploidy_by_sample = ch_input.map { meta, reads -> tuple(meta.sample, (meta.ploidy ?: (meta.n_hap ?: 2))) }
    ch_hgsize_by_sample = ch_input.map { meta, reads -> tuple(meta.sample, meta.haploid_genome_size) }   // column>param, or null
    ch_input            = ch_input.map { meta, reads -> tuple(meta.findAll { k, v -> !(k in ['ploidy', 'haploid_genome_size']) }, reads) }
    ch_versions = Channel.empty()
    ch_gxdb_dir = Channel.empty()

    // ── Ensure the NCBI taxdump once, up front. Idempotent: the module skips the download
    //    when names.dmp/nodes.dmp are already present at the target path. Needed by
    //    RESOLVE_TAXONOMY and (when enabled) the decon evidence branch.
    DOWNLOAD_TAXDUMP(
        Channel.value( file(params.diamond_taxdump_dir) ),
        Channel.value( (params.diamond_force ?: false) as boolean )
    )
    ch_taxdump     = DOWNLOAD_TAXDUMP.out.taxdump_dir
    ch_taxdump_dir = ch_taxdump          // reused by the decon evidence branch

    // ── 4b-i: resolve organism taxonomy (taxid -> name / kingdom / BUSCO lineage) ──
    RESOLVE_TAXONOMY(
        ch_input
            .map    { meta, reads -> meta.taxid }
            .filter { it != null }
            .unique()
            .combine( ch_taxdump )          // pair each distinct taxid with the taxdump
    )

    ch_taxonomy = RESOLVE_TAXONOMY.out.tsv
        .map { taxid, tsv -> tsv }
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            tuple(row.taxid.toString(),
                  [ name          : organismName(row),
                    kingdom       : kingdomFlag(row),
                    busco_lineage : buscoLineageFor(row),
                    genetic_code  : geneticCodeFor(row),
                    telomere_motif: telomereMotifFor(row) ]) }

    // Per-sample identity side-channel: (sample, {taxid, name, kingdom, busco_lineage}).
    ch_sample_identity = ch_input
        .map { meta, reads -> tuple(meta.taxid?.toString(), meta.sample, meta.species) }
        .filter { taxid, sample, override -> taxid != null }
        .combine(ch_taxonomy, by: 0)                       // many samples per taxid
        .map { taxid, sample, override, tax ->
            tuple(sample, [ taxid         : taxid,
                            name          : override ?: tax.name,
                            kingdom       : tax.kingdom,
                            busco_lineage : tax.busco_lineage,
                            genetic_code  : (params.mitochondria_genetic_code ?: tax.genetic_code),
                            telomere_motif: (params.telomere_motif        ?: tax.telomere_motif) ]) }

    ch_sample_identity.subscribe { sample, tax ->
        log.info "  resolved '${sample}': taxid=${tax.taxid}  name='${tax.name}'  kingdom=${tax.kingdom}  busco=${tax.busco_lineage}  gcode=${tax.genetic_code}  telo=${tax.telomere_motif}"
    }

    // Per-taxid resolved mito genetic code (param override wins over taxonomy-derived).
    ch_gcode_by_taxid = ch_taxonomy
        .map { taxid, tax -> tuple(taxid, (params.mitochondria_genetic_code ?: tax.genetic_code)) }

    // Per-taxid resolved telomere motif (param override wins over taxonomy-derived).
    ch_telo_by_taxid = ch_taxonomy
        .map { taxid, tax -> tuple(taxid, (params.telomere_motif ?: tax.telomere_motif)) }

    // Per-taxid organelle-assembly specs for GetOrganelle (short-read path). Types from the
    // resolved kingdom (plant => plastid + mito; else mito); per-organelle -R/-k/coverage
    // defaults from the taxonomy helpers, each param-overridable. Side-channel (not meta), so
    // tuning these re-runs only GetOrganelle.
    ch_organelle_by_taxid = ch_taxonomy
        .map { taxid, tax ->
            def types = params.getorganelle_organelle_types
                            ? params.getorganelle_organelle_types.tokenize(',')*.trim()
                            : organelleTypesFor(tax)
            def specs = types.collect { t ->
                [ type      : t,
                  recursion : (params.getorganelle_recursion ?: getorganelleRecursionFor(t)),
                  kmers     : (params.getorganelle_kmers     ?: getorganelleKmersFor(t)),
                  coverage  : (params.getorganelle_coverage  ?: getorganelleCoverageFor(t)),
                  word_size : params.getorganelle_word_size ]
            }
            tuple(taxid, specs)
        }

    // Bundled per-sample hifiasm inputs: telomere motif (from resolved identity), organism ploidy,
    // and the haploid-size override. One join-by-sample instead of three separate side inputs.
    ch_hifiasm_traits = ch_sample_identity
        .map { sample, id -> tuple(sample, id.telomere_motif) }
        .join( ch_ploidy_by_sample )
        .join( ch_hgsize_by_sample )
        .map { sample, telo, ploidy, hg -> tuple(sample, [ telomere_motif: telo, ploidy: ploidy, haploid_genome_size: hg ]) }

    /*
    ========================================================================================
        STEP 0: Setup Decontamination Databases (if enabled) & BUSCO Database
        Runs in parallel with BAM conversion and assembly
        Only executes if decontamination is requested
    ========================================================================================
    */
    if (params.run_decon_contigs || (capabilities.scaffold && params.run_decon_scaffolds)) {
        SETUP_DECONTAM_DBS(ch_taxdump)
        
        // Store outputs for later use
        ch_gxdb_dir = SETUP_DECONTAM_DBS.out.gxdb_dir
        ch_diamond_db = SETUP_DECONTAM_DBS.out.diamond_db
    }
    
    // BUSCO databases — one download per DISTINCT lineage across all samples
    // (storeDir dedupes; the per-lineage tasks run in parallel).
    ch_busco_lineages = ch_taxonomy
        .filter { qc_on }
        .map { taxid, tax -> tax.busco_lineage }
        .unique()
    ch_busco_downloads = Channel.empty()
    if (qc_on) {
        DOWNLOAD_BUSCO_DB(ch_busco_lineages)
        ch_busco_downloads = DOWNLOAD_BUSCO_DB.out.db
    }

    // ch_busco_db is now a VALUE-channel MAP:  taxid -> busco_lineage (a String).
    //   * value channel  -> broadcasts unchanged to all 13 ASSEMBLY_QC calls
    //   * the .combine on DOWNLOAD_BUSCO_DB.out forces each lineage to finish
    //     downloading before BUSCO reads it from the shared --download_path
    //   * .reduce collapses to ONE map, so the channel emits once (after all
    //     downloads) and broadcasts; multiplicity is what would otherwise break
    //     the fan-out to 13 subworkflow calls.
    // NOTE: many taxids can share one lineage -> combine(by:0) (one-to-many), NOT join.
    ch_busco_db = ch_taxonomy
        .map { taxid, tax -> tuple(tax.busco_lineage, taxid) }
        .combine( ch_busco_downloads.map { db -> tuple(db.name, db) }, by: 0 )
        .map { lineage, taxid, db -> [ (taxid): lineage ] }
        .reduce([:]) { acc, m -> acc + m }
        .map { m -> m.sort { a, b -> a.key <=> b.key } }        // deterministic key order -> stable cache hash

    /*
    ========================================================================================
        STEP 0b: Find Closest Mitochondrial Reference (NCBI)
        Runs at pipeline start — no dependencies on reads or assembly
    ========================================================================================
    */
    // Mito reference per DISTINCT resolved species name, among HiFi samples only
    // (short-read samples use mitofinder / the organelle `other` branch — no MitoHiFi ref).
    // Runs after taxonomy resolution now (was a global no-dependency call).
    ch_mito_ref_todo = ch_input
        .filter { meta, reads -> meta.hifi }
        .map    { meta, reads -> tuple(meta.taxid?.toString(), meta.sample) }
        .filter { taxid, sample -> taxid != null }
        .combine( ch_taxonomy, by: 0 )                 // (taxid, sample, tax) — many samples per taxid
        .map    { taxid, sample, tax -> tuple(taxid, tax.name) }
        .filter { taxid, name -> name != null }        // no resolved name -> no reference
        .unique()                                      // distinct (taxid, name)
    ch_mito_ref_by_taxid = Channel.empty()
    if (capabilities.hifi) {
        FIND_MITO_REFERENCE(ch_mito_ref_todo)
        // Per-taxid reference for the organelle step: (taxid, ref_fasta, ref_gb).
        ch_mito_ref_by_taxid = FIND_MITO_REFERENCE.out.ref_fasta
            .join(FIND_MITO_REFERENCE.out.ref_gb)
    }

    /*
    ========================================================================================
        STEPS 1-3: Prepare reads and estimate genome size
    ========================================================================================
    */
    
    READ_PREPARATION(ch_input, ch_ploidy_by_sample, capabilities)
    ch_reads_all = READ_PREPARATION.out.reads
    ch_qc_reads = READ_PREPARATION.out.qc_reads
    ch_shortread_reads = READ_PREPARATION.out.shortread
    ch_versions = ch_versions.mix(READ_PREPARATION.out.versions)


    /*
    ========================================================================================
        STEP 3b: Assemble Mitochondrial Genome (MitoHiFi)
        Runs concurrently with HIFIASM — both consume HiFi reads independently
        Depends on: BAM_TO_FASTQ + FIND_MITO_REFERENCE
    ========================================================================================
    */
    ORGANELLE(
        ch_reads_all.map { meta, hifi_fastq, hic1, hic2, sr1, sr2 -> tuple(meta, hifi_fastq, sr1, sr2) },
        ch_mito_ref_by_taxid,
        ch_gcode_by_taxid,
        ch_organelle_by_taxid, capabilities
    )
    ch_versions = ch_versions.mix(ORGANELLE.out.versions)

    /*
    ========================================================================================
        STEP 4: Assemble contigs — assembler selector (hifiasm | spades)
        hifiasm: HiFi(+Hi-C) -> hap1+hap2 (diploid) / primary (haploid)
        spades:  PE short reads -> one collapsed 'primary' assembly
    ========================================================================================
    */

    CONTIG_ASSEMBLY(
        ch_reads_all,
        ch_hifiasm_traits,
        READ_PREPARATION.out.genome_size.map { meta, f -> tuple(meta.sample, f) }, capabilities
    )
    ch_versions = ch_versions.mix(CONTIG_ASSEMBLY.out.versions)

    // Fork sample-level assembly into per-hap contigs: diploid -> hap1+hap2, haploid/spades -> one primary.
    CONTIG_ASSEMBLY.out.assemblies
        .flatMap { meta, fastas ->
            def hmetas = forkHaplotypeMeta(meta)           // [hap1, hap2]  or  [primary]
            def fs = (fastas instanceof List) ? fastas.sort { it.name } : [fastas]
            if (fs.size() != hmetas.size()) error "${meta.sample}: expected ${hmetas.size()} assembly files, received ${fs.size()}"
            [hmetas, fs].transpose().collect { hm, fa -> tuple(hm, fa) }
        }
        .set { ch_contigs }                                // per-haplotype tuple(meta, fasta)

    /*
    ====================================================================================
        STEP 4a: Remove Organelle Contigs from Nuclear Assemblies (mito + plastid)
        Runs on ALL contigs (HiFi and short-read), each baited by that sample's own
        organelle assemblies (ORGANELLE.out.assemblies). Samples that produced no organelle
        pass through untouched. The filtered contigs are then forked by read type into the
        HiFi and short-read conditioning paths.
    ====================================================================================
    */
    // Per-sample bait = ALL organelle assemblies for that sample (mito + plastid), emitted
    // one-per-sample by ORGANELLE with branch-local gathers (see workflows/organelle.nf) so
    // one sample can't stall the rest. An empty bait list means the sample produced no
    // organelle: those bypass FILTER_ORGANELLE untouched, so a failed organelle assembly
    // still can't drop a whole nuclear assembly.
    CONTIG_REFINEMENT(ch_contigs, ORGANELLE.out.baits, READ_PREPARATION.out.hifi,
        ch_shortread_reads, ch_gxdb_dir, capabilities)
    ch_decontaminated_contigs = CONTIG_REFINEMENT.out.assembly
    ch_organelle_filtered = CONTIG_REFINEMENT.out.organelle_filtered
    ch_shortread_conditioned = CONTIG_REFINEMENT.out.shortread_conditioned
    ch_versions = ch_versions.mix(CONTIG_REFINEMENT.out.versions)
    HIC_SCAFFOLDING(ch_decontaminated_contigs, READ_PREPARATION.out.hifi,
        READ_PREPARATION.out.hic, ch_telo_by_taxid, ch_gxdb_dir, capabilities)
    ch_final_assembly = HIC_SCAFFOLDING.out.assembly
    ch_shortread_finished = HIC_SCAFFOLDING.out.shortread
    ch_final_scaffolds_round2 = HIC_SCAFFOLDING.out.round2
    ch_scaffold_round2_agp = HIC_SCAFFOLDING.out.round2_agp
    ch_all_bam_metrics = HIC_SCAFFOLDING.out.bam_metrics
    ch_all_pairs_metrics = HIC_SCAFFOLDING.out.pairs_metrics
    ch_teloclip_stats_for_report = HIC_SCAFFOLDING.out.teloclip_stats
    ch_versions = ch_versions.mix(HIC_SCAFFOLDING.out.versions)
    // =========================================================================
    //  FINALIZE ASSEMBLY — now uses ch_final_assembly (post-teloclip if enabled)
    // =========================================================================
    // Harmonize scaffold names across same-species long-read assemblies (>= 2) before
    // finalizing, so homologous chromosomes share names and FINAL_VIZ inherits them.
    // Short-read assemblies are never harmonized (harmonize_scaffolds.nf branches them out),
    // and mixing them into the INPUT here would hold its cross-sample groupTuple open until
    // the short-read path finished -- blocking FINALIZE and everything after it for every
    // long-read assembly. They bypass it and rejoin at the output with the same sentinel.
    HARMONIZE_SCAFFOLDS(ch_final_assembly, ch_harmonize_script, capabilities)
    ch_versions = ch_versions.mix(HARMONIZE_SCAFFOLDS.out.versions)

    // ---- break chimeric scaffolds, if asked -------------------------------------------
    // Between harmonization and finalization, and that is the only possible place:
    // harmonization emits a name MAP and FINALIZE_ASSEMBLY applies it, so this is the last
    // point at which the FASTA is in original coordinates and the composite is one record,
    // and the first at which the concordance vote exists to justify cutting it.
    //
    // ---- chimeric scaffold detection, confirmation and repair ------------------------
    // Lifted into workflows/chimera.nf: main.nf exceeded Groovy's 65,535-character
    // compiled-unit limit. Workflow nesting is part of task identity.
    //
    // Detection runs whenever harmonization did, independently of chimera_break -- the
    // candidates and called joins are the evidence a cut is justified by, and they are worth
    // having on a run that cuts nothing. Breaking is gated separately inside.
    CHIMERA(
        HARMONIZE_SCAFFOLDS.out.assemblies,
        ch_shortread_finished,
        HIC_SCAFFOLDING.out.round1_agp,
        ch_scaffold_round2_agp,
        HARMONIZE_SCAFFOLDS.out.ref_pafs_by_id,
        HARMONIZE_SCAFFOLDS.out.chimera_candidates,
        HARMONIZE_SCAFFOLDS.out.ref_name_map,
        HIC_SCAFFOLDING.out.contig_pairs,
        ch_telo_by_taxid, capabilities )
    ch_versions = ch_versions.mix(CHIMERA.out.versions)

    // ch_pre_finalize is emitted rather than assigned here: its default carries the
    // short-read-only branch (no harmonization name map), and separating that default from
    // the BREAK_CHIMERAS override is what previously dropped those assemblies silently.
    ch_pre_finalize   = CHIMERA.out.pre_finalize
    ch_chimeric_joins = CHIMERA.out.called

    ch_name_map_files = HARMONIZE_SCAFFOLDS.out.assemblies
        .map { meta, fa, nm -> nm }
        .filter { nm -> !nm.name.startsWith('NO_') }
        .collect()
    COLLECT_NAME_MAPS( ch_name_map_files )
    ch_name_map_for_report = COLLECT_NAME_MAPS.out.map.ifEmpty( file('NO_NAMEMAP') )

    FINALIZE_ASSEMBLY(ch_pre_finalize)
    ch_finalized_assembly = FINALIZE_ASSEMBLY.out.assembly
    ASSEMBLY_RUN_SUMMARY(ch_input_status,
        ch_finalized_assembly.map { meta, fa -> [sample: meta.sample, id: meta.id, name: fa.name] }.collect(flat: false))

    // =========================================================================
    //  PanGenome Assembly - combine all same species chromosome level assemblies into a pangenome
    // =========================================================================

    // Pangenome graph (minigraph-cactus) per species, from the gated finalized assemblies.
    ch_finalized_with_fai = FINALIZE_ASSEMBLY.out.assembly.join(FINALIZE_ASSEMBLY.out.fai)

    /*
    ========================================================================================
        BUILD MERYL DATABASE (Once per sample)
        Runs in parallel with trimming and assembly
        Reused across ALL assembly QC steps for dramatic speedup
    ========================================================================================
    */
    ch_meryl_db = Channel.empty()
    if (qc_on || params.run_pangenome) {
        BUILD_MERYL_DB(ch_qc_reads)
        ch_meryl_db = BUILD_MERYL_DB.out.meryl_db
        ch_versions = ch_versions.mix(BUILD_MERYL_DB.out.versions)
    }

    // resolved organism name per taxid (RESOLVE_TAXONOMY -> ch_taxonomy; tax.name is the
    // taxid-derived species, e.g. 373251 -> "Spratelloides delicatulus")
    ch_species_by_taxid = ch_taxonomy.map { taxid, tax -> tuple(taxid.toString(), tax.name) }

    PANGENOME(
        ch_finalized_with_fai,
        HARMONIZE_SCAFFOLDS.out.reference_id,
        ch_species_by_taxid,
        HARMONIZE_SCAFFOLDS.out.report_by_taxid,
        ch_meryl_db
    )
    ch_versions = ch_versions.mix(PANGENOME.out.versions)

    ch_pangenome_report_for_report = PANGENOME.out.report
        .map { taxid, md -> md }
        .ifEmpty( file('NO_PANGENOME') )

    /*
    ========================================================================================
        STEP 14: Finalization
    ========================================================================================
    */
    // 1. Contact Maps for Final Assemblies
    //    REPLACE: HIC_SCAFFOLDING.out.filled → ch_final_assembly
    ch_final_contact_maps = Channel.empty()
    if (post_on && capabilities.hic && params.run_final_contact_maps) {
        FINAL_HIC_MAPS(
            ch_finalized_assembly,
            READ_PREPARATION.out.hic,
            ch_compartments_script,
            ch_tad_book_script
        )
        ch_all_bam_metrics   = ch_all_bam_metrics.mix(FINAL_HIC_MAPS.out.bam_metrics)
        ch_all_pairs_metrics = ch_all_pairs_metrics.mix(FINAL_HIC_MAPS.out.pairs_metrics)
        ch_final_contact_maps = FINAL_HIC_MAPS.out.contact_maps
    }

    // String-id view of the finalized per-hap assemblies for the leaf viz/QUAST steps.
    // These only need an id (for pairing, labels, output naming) — meta itself isn't used.
    ch_finalized_assembly
        .map { meta, fa -> tuple(meta.id, fa) }
        .set { ch_final_by_id }

    ch_finalized_assembly
        .map { meta, fa -> tuple(meta.taxid?.toString(), meta.id, fa) }
        .combine( ch_telo_by_taxid, by: 0 )
        .map { taxid, id, fa, telo -> tuple(id, fa, telo) }
        .set { ch_final_by_id_telo }

    // 3. Dotplots of final assemblies vs each other (hap1 vs hap2)
    /*
    ========================================================================================
        PAIRWISE GENOME ALIGNMENTS
        Generates all pairwise alignments for dotplot visualization grid
    ========================================================================================
    */
    // Defaulted BEFORE the branch: a variable assigned only inside an if-block is not
    // visible after it in a workflow body. That is the "No such variable" failure the
    // chimera wiring hit, and these two are read by REPORTING further down.
    ch_pairwise_summary    = Channel.empty()
    ch_telomere_for_report = Channel.empty()
    ch_viz_dotplot         = Channel.empty()
    ch_viz_riparian        = Channel.empty()
    ch_viz_tidk            = Channel.empty()

    if (post_on) {
        FINAL_VIZ(ch_final_by_id, ch_final_by_id_telo, ch_dotplot_script, ch_riparian_script)
        ch_versions            = ch_versions.mix(FINAL_VIZ.out.versions)
        ch_pairwise_summary    = FINAL_VIZ.out.pairwise_summary
        ch_telomere_for_report = FINAL_VIZ.out.telomere_summary
        ch_viz_dotplot         = FINAL_VIZ.out.dotplot
        ch_viz_riparian        = FINAL_VIZ.out.riparian
        ch_viz_tidk            = FINAL_VIZ.out.tidk_plot
    }
    else {
        log.info "[INFO] run_post_assembly = false: skipping FINAL_VIZ"
    }

    // 6. NCBI output files for GenBank submission (if enabled)
    
    /*
    ========================================================================================
        QC Steps
    ========================================================================================
    */

    /*
    ========================================================================================
        Sequencing QC
    ========================================================================================
    */
    /*
    ========================================================================================
        QC Raw Hi-C Reads
    ========================================================================================
    */
    // Register only read-QC branches supported by accepted inputs.
    if (qc_on && capabilities.hic) {
        HIC_QC_RAW(READ_PREPARATION.out.hic_raw, "raw")
        HIC_QC_TRIMMED(READ_PREPARATION.out.hic_trimmed, "trimmed")
        ch_versions = ch_versions.mix(HIC_QC_RAW.out.versions)
    }
    if (qc_on && capabilities.hifi) {
        HIFI_QC(READ_PREPARATION.out.hifi)
        ch_versions = ch_versions.mix(HIFI_QC.out.versions)
    }
    if (qc_on && capabilities.shortread) {
        SHORTREAD_QC_RAW(ch_input.filter { meta, reads -> meta.shortread }
            .map { meta, reads -> tuple(meta, reads.sr_r1, reads.sr_r2) }, "raw")
        ch_versions = ch_versions.mix(SHORTREAD_QC_RAW.out.versions)
        if (params.run_shortread_trim) {
            SHORTREAD_QC_TRIMMED(READ_PREPARATION.out.shortread, "trimmed")
        }
    }
    /*
    ========================================================================================
        Assembly QC
    ========================================================================================
    */

    // Stage every checkpoint's per-hap assembly into one labeled channel for QC_PHASE.
    // Guards below decide which stages are QC'd (lifted from the old per-call ifs); 'final' is unconditional.
    ch_staged_assemblies = Channel.empty()
    if (run_all_qc) ch_staged_assemblies = ch_staged_assemblies.mix( ch_contigs.map { m, f -> tuple(m, 'initial', f) } )
    if (run_all_qc) ch_staged_assemblies = ch_staged_assemblies.mix( ch_organelle_filtered.map { m, f -> tuple(m, 'organelle_filtered', f) } )
    if (run_all_qc) ch_staged_assemblies = ch_staged_assemblies.mix( CONTIG_REFINEMENT.out.purged.map { m, f -> tuple(m, 'purged', f) } )
    if (run_all_qc) ch_staged_assemblies = ch_staged_assemblies.mix( ch_shortread_conditioned.map { m, f -> tuple(m, 'redundans', f) } )
    if (run_all_qc && params.run_inspector_contigs) ch_staged_assemblies = ch_staged_assemblies.mix( CONTIG_REFINEMENT.out.corrected.map { m, f -> tuple(m, 'contig_corrected', f) } )
    if (run_all_qc && params.run_decon_contigs) ch_staged_assemblies = ch_staged_assemblies.mix( CONTIG_REFINEMENT.out.decontaminated.map { m, f -> tuple(m, 'contig_decontam', f) } )
    if (run_all_qc) ch_staged_assemblies = ch_staged_assemblies.mix( HIC_SCAFFOLDING.out.round1.map { m, f -> tuple(m, 'scaffold', f) } )
    if (run_all_qc && params.run_inspector_scaffolds) ch_staged_assemblies = ch_staged_assemblies.mix( HIC_SCAFFOLDING.out.corrected.map { m, f -> tuple(m, 'scaffold_corrected', f) } )
    if (run_all_qc && params.run_decon_scaffolds) ch_staged_assemblies = ch_staged_assemblies.mix( HIC_SCAFFOLDING.out.decontaminated.map { m, f -> tuple(m, 'scaffold_decontam', f) } )
    if (run_all_qc && params.run_scaffold_round2) ch_staged_assemblies = ch_staged_assemblies.mix( ch_final_scaffolds_round2.map { m, f -> tuple(m, 'scaffold_round2', f) } )
    if (run_all_qc) ch_staged_assemblies = ch_staged_assemblies.mix( HIC_SCAFFOLDING.out.filled.map { m, f -> tuple(m, 'gap_filled', f) } )
    if (run_all_qc && params.run_teloclip_extend) ch_staged_assemblies = ch_staged_assemblies.mix( HIC_SCAFFOLDING.out.extended.map { m, f -> tuple(m, 'teloclip', f) } )
    if (run_all_qc && chimera_on) ch_staged_assemblies = ch_staged_assemblies.mix( CHIMERA.out.broken.map { m, f, nm -> tuple(m, 'chimera_broken', f) } )
    if (qc_on) ch_staged_assemblies = ch_staged_assemblies.mix( ch_finalized_assembly.map { m, f -> tuple(m, 'final', f) } )

    // Everything derived from QC_PHASE goes with it: ch_final_busco, ch_final_busco_table,
    // SNAIL_PLOT_FINAL (which joins the BUSCO full table) and REPORTING's metrics/plots.
    ch_final_busco       = Channel.empty()
    ch_final_busco_table = Channel.empty()
    ch_qc_metrics        = Channel.empty()
    ch_qc_plots          = Channel.empty()
    ch_qc_report_html    = Channel.empty()

    if (qc_on) {
    QC_PHASE(
        ch_staged_assemblies,
        ch_qc_reads,
        ch_meryl_db,
        ch_busco_db,
        ch_all_bam_metrics,
        ch_all_pairs_metrics,
        ch_compile_qc_script,
        ch_assembly_report_script, capabilities
    )
    ch_versions = ch_versions.mix(QC_PHASE.out.versions)
    }
    else {
        log.info "[INFO] qc_mode = 'none': skipping QC_PHASE and everything derived from it"
    }

    /*
    ========================================================================================
        Generate Optional Decontamination Evidence
    ========================================================================================
    */
    // BlobTools evidence remains unsupported pending its separate rebuild.

    // Guarded where they stand rather than by extending the QC_PHASE block: the span
    // between them contains the scaffold decontamination report branch, which is not a QC
    // artifact and must not become conditional on qc_mode. The Channel.empty() defaults
    // above stand when QC is off.
    if (qc_on) {
        ch_final_busco       = QC_PHASE.out.final_busco
        ch_qc_metrics        = QC_PHASE.out.metrics
        ch_qc_plots          = QC_PHASE.out.plots
        ch_qc_report_html    = QC_PHASE.out.report_html
    }

    /*
    ========================================================================================
    SNAIL PLOTS FOR FINAL ASSEMBLIES
    ========================================================================================
    */
    // Join gap-filled assemblies with their BUSCO results
    // BUSCO output from ASSEMBLY_QC_TELOCLIP is per-haplotype
    if (qc_on) ch_final_busco_table = QC_PHASE.out.final_busco_table

    ch_finalized_assembly
        .join(ch_final_busco_table)
        .map { meta, assembly, full_table ->
            tuple(meta.id, assembly, full_table, "final")
        }
        .set { ch_snail_plot_final_input }

    // QC-dependent: it joins ch_final_busco_table, which only exists when QC_PHASE ran.
    // Its output is read by REPORTING, so it needs the same defaulted-channel treatment as
    // the other gated processes: reading a gated process output directly is what produced
    // the "Access to .out is undefined" failure three times in a row here.
    ch_snail_final = Channel.empty()
    if (qc_on) {
        SNAIL_PLOT_FINAL(ch_snail_plot_final_input)
        ch_snail_final = SNAIL_PLOT_FINAL.out.snail
    }
    /*
    ========================================================================================
        COVERAGE BOOK - HiFi coverage visualization for final assemblies
    ========================================================================================
    */

    // Combine gap-filled assemblies with HiFi reads for coverage book
    ch_finalized_assembly
        .map { meta, final_fa -> tuple(meta.sample, meta, final_fa) }
        .combine(
            READ_PREPARATION.out.hifi.map { meta, hifi_fastq -> tuple(meta.sample, hifi_fastq) },
            by: 0
        )
        .map { sample, meta, final_fa, hifi_fastq ->
            tuple(meta, final_fa, hifi_fastq)
        }
        .set { ch_coverage_book_input }

    if (post_on) COVERAGE_BOOK(ch_coverage_book_input, ch_coverage_book_script)

    // =========================================================================
    //  SUMMARY REPORT  (manifest assembly + report -> workflows/reporting.nf)
    // =========================================================================
    // ---- chimera tables for the report ---------------------------------------------
    // Built inside CHIMERA, which owns the channels they come from.
    ch_chimera_candidates_rpt = CHIMERA.out.candidates_for_rpt
    ch_chimera_joins_rpt      = CHIMERA.out.joins_for_rpt
    ch_chimera_evidence_rpt   = CHIMERA.out.evidence_for_rpt
    ch_chimera_figures_rpt    = CHIMERA.out.figures_for_rpt

    // Consumes outputs of both QC_PHASE and FINAL_VIZ, so it needs both.
    if (qc_on && post_on) {
    REPORTING(
        ch_finalized_assembly,
        ch_snail_final,
        ch_final_contact_maps,
        ch_viz_dotplot,
        ch_viz_riparian,
        ch_viz_tidk,
        ch_qc_metrics,
        ch_qc_plots,
        ch_qc_report_html,
        ORGANELLE.out.annotation,
        ORGANELLE.out.stats,
        ORGANELLE.out.circular_map,
        READ_PREPARATION.out.genome_results,
        READ_PREPARATION.out.genome_size,
        ch_sample_identity,
        ch_input,
        ch_ploidy_by_sample,
        ch_telomere_for_report,
        ch_pairwise_summary,
        ch_teloclip_stats_for_report,
        ch_pangenome_report_for_report,
        ch_name_map_for_report,
        ch_chimera_candidates_rpt,
        ch_chimera_joins_rpt,
        ch_chimera_evidence_rpt,
        ch_chimera_figures_rpt,
        ch_versions,
        ch_summary_report_script
    )
    }
    else {
        log.info "[INFO] skipping REPORTING: it needs both QC and the post-assembly outputs"
    }
}
/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """\
        Pipeline completed!
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Results   : ${params.outdir}
        """
        .stripIndent()
}
