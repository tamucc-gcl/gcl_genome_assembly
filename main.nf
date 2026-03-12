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
params.sample_sheet = null
params.outdir = './results'
params.publish_dir_mode = 'link'

// Assembly parameters
////   Overlap/Error correction:
params.hifiasm_k = 51 // k-mer length (must be <64) [51]
params.hifiasm_w = 51 // minimizer window size [51]
params.hifiasm_f = 37 // number of bits for bloom filter; 0 to disable [37]
params.hifiasm_D = 5.0 //Drop k-mers occurring >FLOAT*coverage times [5.0]
params.hifiasm_N = 100 //consider up to max(-D*coverage,-N) overlaps for each oriented read [100]
params.hifiasm_r = 3 //round of correction [3]
params.hifiasm_z = 0 //length of adapters that should be removed [0]
params.hifiasm_maxKOCC = 2000 //employ k-mers occurring <INT times to rescue repetitive overlaps [2000]
params.hifiasm_hgSize = 'auto' //estimated haploid genome size used for inferring read coverage [auto] (INT(k, m or g))
//// Assembly
params.hifiasm_a = 4 // round of assembly cleaning [4]
params.hifiasm_m = 10000000 // pop bubbles of <INT in size in contig graphs [10000000]
params.hifiasm_p = 0 // pop bubbles of <INT in size in unitig graphs [0]
params.hifiasm_n = 3 // remove tip unitigs composed of <=INT reads [3]
params.hifiasm_x = 0.8 // max overlap drop ratio [0.8]
params.hifiasm_y = 0.2 // min overlap drop ratio [0.2]
params.hifiasm_u = 1 // post-join step for contigs which may improve N50; 0 to disable; 1 to enable
params.hifiasm_homCov = 'auto' // homozygous read coverage [auto]  Int
params.hifiasm_lowQ = 70 // output contig regions with >=INT% inconsistency in BED format; 0 to disable [70]
params.hifiasm_bCov =  0 // break contigs at positions with <INT-fold coverage; work with '--m-rate'; 0 to disable [0]
params.hifiasm_hCov = -1 // break contigs at positions with >INT-fold coverage; work with '--m-rate'; -1 to disable [-1]
params.hifiasm_mRate = 0.75 //break contigs at positions with >INT-fold coverage; work with '--m-rate'; -1 to disable [-1] only work with '--b-cov' or '--h-cov'[0.75]
params.hifiasm_primary = false // output a primary assembly and an alternate assembly
params.hifiasm_ctgN = 3 // remove tip contigs composed of <=INT reads [3]
////Purge-dups:
params.hifiasm_l = 3 // purge level. 0: no purging; 1: light; 2/3: aggressive [0 for trio; 3 for unzip]
params.hifiasm_s = 0.55 // similarity threshold for duplicate haplotigs in read-level [0.75 for -l1/-l2, 0.55 for -l3]
params.hifiasm_O = 1 // min number of overlapped reads for duplicate haplotigs [1]
params.hifiasm_purgeMax = 'auto' // coverage upper bound of Purge-dups [auto] Int
params.hifiasm_nHaplotypes = 2 // number of haplotypes to identify - only works with 1/2
////Hi-C-partition:
params.hifiasm_useHiC = true //use hi-c reads in contig assembly?
params.hifiasm_sBase = 0.5 //similarity threshold for homology detection in base-level; -1 to disable [0.5]; -s for read-level (see <Purge-dups>)
params.hifiasm_nWeight = 3 //rounds of reweighting Hi-C links [3]
params.hifiasm_nPerturb = 10000 // rounds of perturbation [10000]
params.hifiasm_fPerturb = 0.1 // fraction to flip for perturbation [0.1]
params.hifiasm_lMSjoin = 500000 // detect misjoined unitigs of >=INT in size; 0 to disable [500000]
////Dual-Scaffolding:
params.hifiasm_dualScaf = false // output scaffolding
params.hifiasm_scafGap = 3000000 // max gap size for scaffolding [3000000]
////Telomere-identification:
// see telomere specific parameters below (used across modules)
params.hifiasm_teloP = 1 // non-telomeric penalty [1]
params.hifiasm_teloD = 2000 // max drop [2000]
params.hifiasm_teloS = 500 // min score for telomere reads [500]

// Mitochondrial genome assembly (MitoHiFi)
params.mitohifi_species        = null        // REQUIRED: Scientific name (e.g., "Spratelloides delicatulus")
params.mitohifi_genetic_code   = 2           // NCBI genetic code: 2=vertebrate mito, 5=invertebrate mito
params.mitohifi_ref_min_length = 14000       // Min reference mitogenome length for findMitoReference.py
params.mitohifi_perc_identity  = 50          // Min percent identity for read filtering
params.mitohifi_cov_cutoff     = 'auto'      // Coverage cutoff (auto or integer)
params.mitohifi_bloom_filter   = false       // Use bloom filter for large read sets
params.mitohifi_circular_map_dpi = 300

// Mito contig filtering (applied to HIFIASM output)
params.mitohifi_filter_min_identity = 90     // Min alignment identity to call a contig mitochondrial
params.mitohifi_filter_min_coverage = 50     // Min query coverage to call a contig mitochondrial


// Purge duplicates parameters
params.run_purge_dups = false  // Default off - enable with --run_purge_dups true

// Assembly QC parameters
params.busco_lineage = 'actinopterygii_odb10'
params.busco_downloads = '/work/birdlab/databases/busco'
params.merqury_k = 21  // k-mer size for meryl database

// Misassembly correction parameters (Inspector) - CONTIGS
params.inspector_run_on_contigs = true
params.inspector_contig_skip_baseerror = true
params.inspector_contig_min_depth = null  // default: 20% of average depth
params.inspector_contig_min_contig_length = 10000
params.inspector_contig_min_contig_length_assemblyerror = 1000000
params.inspector_contig_min_assembly_error_size = 50
params.inspector_contig_max_assembly_error_size = 4000000

// Misassembly correction parameters (Inspector) - SCAFFOLDS
params.inspector_run_on_scaffolds = true
params.inspector_scaffold_skip_baseerror = true
params.inspector_scaffold_min_depth = null  // default: 20% of average depth
params.inspector_scaffold_min_contig_length = 10000
params.inspector_scaffold_min_contig_length_assemblyerror = 10000000  // Higher threshold for scaffolds
params.inspector_scaffold_min_assembly_error_size = 50
params.inspector_scaffold_max_assembly_error_size = 10000000  // Higher threshold for scaffolds

// Hi-C mapping and QC parameters
params.bwa_mem2_hic_args = "" // Additional arguments for bwa-mem2 Hi-C mapping
params.hic_coverage_window = 100000
params.hic_min_mapq = 30
params.hic_resolutions = "2500000,1000000,500000,250000,100000,50000,10000"
params.hic_base_bin = "10000"
params.hic_plot_resolutions = "1000000,500000,250000,100000"
params.hic_balance = true
params.hic_min_mapq_raw = 30
params.hic_min_mapq_filtered = 1
params.scaffold_min_size = 0  // 0 = include all scaffolds in contact maps; set to e.g. 100000 to include only scaffolds >100kb

// Hi-C scaffolding parameters (YaHS) - ROUND 1
params.yahs_round1_min_contig_length = 20000
params.yahs_round1_min_mapq = 10
params.yahs_round1_resolutions = '20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000'
params.yahs_round1_rounds_per_resolution = null
params.yahs_round1_enzyme = null
params.yahs_round1_no_contig_ec = false
params.yahs_round1_no_scaffold_ec = false

// Hi-C scaffolding parameters (YaHS) - ROUND 2
params.yahs_round2_min_contig_length = 100000
params.yahs_round2_min_mapq = 20
params.yahs_round2_resolutions = '50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000'
params.yahs_round2_rounds_per_resolution = null
params.yahs_round2_enzyme = null
params.yahs_round2_no_contig_ec = false
params.yahs_round2_no_scaffold_ec = true


// Finalization Parameters
params.make_final_contact_maps = true
params.report_stage = 'final'

// Dotplot Params
params.run_pairwise_alignments = true
params.pairwise_alignment_preset = 'asm5'  // minimap2 preset for assembly vs assembly alignment
params.pairwise_alignment_min_mapq = 5
params.pairwise_alignment_min_aln_bp = 10000
params.pairwise_alignment_mode = 'within_sample'         // 'all' = all pairs, 'within_sample' = hap1 vs hap2 only
params.pairwise_dotplot_width = 10             // Dotplot width in inches
params.pairwise_dotplot_height = 10            // Dotplot height in inches

params.riparian_min_aln_bp = 50000
params.riparian_min_seq_bp = 1000000
params.riparian_alpha      = 0.45
params.riparian_width      = 14
params.riparian_height     = 6

// Compartments parameters - only made when contact maps are made and balanced
params.compartment_resolution = 250000
params.compartment_min_contig_bp = 5000000
params.compartment_max_contigs = 30

// TADs parameters
params.tad_resolution = 50000 // must be in hic_resolutions
params.tad_window_bp = 500000 // must be in hic_resolutions
params.tad_min_contig_bp = 5000000
params.tad_max_contigs = 0    // 0 = no cap

// ============================================================================
// Telomere detection parameters
// ============================================================================
params.telomere_motif = 'CCCTAA'  // telomere motif at 5'-end complement/reverse generated automatically. used by hifiasm (telo-m) & telomere scanning
params.telomere_window = 10000        // Window size (bp) at each sequence end to search
params.telomere_min_repeats = 10      // Minimum consecutive motif repeats required

// Teloclip: telomere extension from soft-clipped HiFi read overhangs
params.run_teloclip_extend     = true       // Master switch for teloclip extension
params.teloclip_min_clip       = 1          // Min overhang bases past contig end
params.teloclip_max_break      = 50         // Max gap between alignment and contig end
params.teloclip_min_anchor     = 100        // Min anchored alignment length
params.teloclip_min_mapq       = 20         // Min mapping quality
params.teloclip_min_overhangs  = 1          // Min supporting overhangs for extension
params.teloclip_max_homopolymer = 500       // Max homopolymer run in overhang

// tidk: telomere identification toolkit (replaces scan_telomeres.py)
params.tidk_explore_minimum    = 5          // Min kmer length for explore
params.tidk_explore_maximum    = 12         // Max kmer length for explore
params.tidk_search_window      = 10000      // Window size for tidk search
params.tidk_plot_height        = 200        // SVG subplot height (px)
params.tidk_plot_width         = 1000       // SVG plot width (px)

// ============================================================================
// genome book parameters
// ============================================================================
params.bin_size = 1000
params.min_len = 1000000
params.min_mapq = 5

// Scaffolding round control
// If not explicitly set, default to true if scaffold correction OR decontamination is enabled
//params.run_scaffold_round2 = null  // Will be computed below if null

// ============================================================================
// Database base directory
// ============================================================================
params.db_base = "/work/birdlab/databases"

// ============================================================================
// FCS-GX database settings (FLATTENED)
// ============================================================================
params.gxdb_dir = "${params.db_base}/fcs-gx"
params.gxdb_profile = 'all'          // 'all' | 'test-only'
params.gxdb_manifest = null          // optional: override manifest URL/path
params.gxdb_force = false            // re-download even if present

// ============================================================================
// DIAMOND / blobtools evidence DB settings (FLATTENED)
// ============================================================================
params.diamond_dmnd = null           // if you have a prebuilt .dmnd, set this
params.diamond_dir = "${params.db_base}/diamond"
params.diamond_name = "proteins"
params.diamond_taxdump_dir = "${params.db_base}/ncbi_taxonomy"
params.diamond_profile = 'custom'    // 'custom' (recommended)
params.diamond_fasta_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
params.diamond_taxonmap_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz"
params.diamond_force = false

// ============================================================================
// Decontamination control (FLATTENED)
// ============================================================================
// When to run decontamination
params.decon_run_on_contigs = true       // Run on initial contig assemblies (before Hi-C mapping)
params.decon_run_on_scaffolds = false     // Run on scaffolded assemblies (after Hi-C scaffolding)

// Core settings
params.decon_source_taxid = 7898          // Actinopterygii; set your species taxid when possible

// Optional: adapter/vector removal
params.decon_run_fcs_adaptor = false      // Requires FCS-adaptor containers; enable only if configured
params.decon_fcsadaptor_mode = 'euk'      // 'euk' or 'prok'
params.decon_container_engine = 'singularity'

// Optional: generate evidence (coverage + taxonomy + blobtools plots)
params.decon_make_blobtools_evidence = true

// ============================================================================
// Evidence generation settings (FLATTENED)
// ============================================================================
params.evidence_map_preset = 'map-hifi'
params.evidence_diamond_max_target_seqs = 1
params.evidence_diamond_evalue = 1e-25
params.evidence_blob_min_contig_len = 1000

// ============================================================================
// BACKWARD COMPATIBILITY LAYER
// ============================================================================
// Reconstruct nested structures for internal workflow code
// This must come AFTER all the flattened parameters are defined above

params.gxdb = [
    dir: params.gxdb_dir,
    profile: params.gxdb_profile,
    manifest: params.gxdb_manifest,
    force: params.gxdb_force
]

params.diamond = [
    dmnd: params.diamond_dmnd,
    dir: params.diamond_dir,
    name: params.diamond_name,
    taxdump_dir: params.diamond_taxdump_dir,
    profile: params.diamond_profile,
    fasta_url: params.diamond_fasta_url,
    taxonmap_url: params.diamond_taxonmap_url,
    force: params.diamond_force
]

params.decon = [
    run_on_contigs: params.decon_run_on_contigs,
    run_on_scaffolds: params.decon_run_on_scaffolds,
    source_taxid: params.decon_source_taxid,
    run_fcs_adaptor: params.decon_run_fcs_adaptor,
    fcsadaptor_mode: params.decon_fcsadaptor_mode,
    container_engine: params.decon_container_engine,
    make_blobtools_evidence: params.decon_make_blobtools_evidence
]

params.evidence = [
    map_preset: params.evidence_map_preset,
    diamond_max_target_seqs: params.evidence_diamond_max_target_seqs,
    diamond_evalue: params.evidence_diamond_evalue,
    blob_min_contig_len: params.evidence_blob_min_contig_len
]

// ============================================================================
// BACKWARD COMPATIBILITY LAYER - YaHS parameters
// ============================================================================
params.yahs_round1 = [
    min_contig_length: params.yahs_round1_min_contig_length,
    min_mapq: params.yahs_round1_min_mapq,
    resolutions: params.yahs_round1_resolutions,
    rounds_per_resolution: params.yahs_round1_rounds_per_resolution,
    enzyme: params.yahs_round1_enzyme,
    no_contig_ec: params.yahs_round1_no_contig_ec,
    no_scaffold_ec: params.yahs_round1_no_scaffold_ec
]

params.yahs_round2 = [
    min_contig_length: params.yahs_round2_min_contig_length,
    min_mapq: params.yahs_round2_min_mapq,
    resolutions: params.yahs_round2_resolutions,
    rounds_per_resolution: params.yahs_round2_rounds_per_resolution,
    enzyme: params.yahs_round2_enzyme,
    no_contig_ec: params.yahs_round2_no_contig_ec,
    no_scaffold_ec: params.yahs_round2_no_scaffold_ec
]

// ============================================================================
// COMPUTE CONDITIONAL PARAMETERS - MUST BE AFTER ALL DEFINITIONS
// ============================================================================
// Determine if second round of scaffolding should run
// Logic: Run round 2 if scaffold correction OR decontamination is enabled (unless explicitly disabled)
if (!params.containsKey('run_scaffold_round2') || params.run_scaffold_round2 == null) {
    // Only compute if not explicitly set by user
    params.run_scaffold_round2 = params.inspector_run_on_scaffolds || params.decon_run_on_scaffolds
}

// Print pipeline header
log.info """\
    =========================================
    GENOME ASSEMBLY PIPELINE
    =========================================
    Sample sheet : ${params.sample_sheet}
    Output dir   : ${params.outdir}
    Inspector    : ${params.inspector_run_on_contigs ? 'Contigs' : ''}${params.inspector_run_on_scaffolds ? ' Scaffolds' : ''}${!params.inspector_run_on_contigs && !params.inspector_run_on_scaffolds ? 'Disabled' : ''}
    Decontamination: ${params.decon.run_on_contigs ? 'Contigs' : ''}${params.decon.run_on_scaffolds ? ' Scaffolds' : ''}${!params.decon.run_on_contigs && !params.decon.run_on_scaffolds ? 'Disabled' : ''}
    Scaffold Round 2: ${params.run_scaffold_round2 ? 'Enabled' : 'Disabled'}
    Gap Filling  : Enabled
    =========================================
    """
    .stripIndent()

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!params.sample_sheet) {
    exit 1, "Sample sheet not specified! Please provide --sample_sheet"
}

/*
========================================================================================
    IMPORT FUNCTIONS
========================================================================================
*/
include { parseSampleSheet } from './functions/parse_sample_sheet.nf'

/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/
include { HIC_QC as HIC_QC_RAW } from './workflows/hic_qc.nf'
include { HIC_QC as HIC_QC_TRIMMED } from './workflows/hic_qc.nf'
include { HIFI_QC } from './workflows/hifi_qc.nf'

// Assembly QC
include { ASSEMBLY_QC as ASSEMBLY_QC_INITIAL } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_MITO_FILTERED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_PURGED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_CONTIG_CORRECTED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_CONTIG_DECONTAM } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_CORRECTED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_DECONTAM } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_ROUND2 } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_GAP_FILLED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_TELOCLIP } from './workflows/assembly_qc.nf'

// HI-C MODULAR WORKFLOWS
include { HIC_QC_FROM_BAM as HIC_QC_FROM_BAM_RAW } from './workflows/hic_qc_from_bam.nf'
include { HIC_QC_FROM_BAM as HIC_QC_FROM_BAM_FILTERED } from './workflows/hic_qc_from_bam.nf'
include { HIC_SCAFFOLD_QC } from './workflows/hic_scaffold_qc.nf'

// DECONTAMINATION MODULAR WORKFLOWS
include { SETUP_DECONTAM_DBS } from './workflows/setup_decontam_dbs.nf'
include { DECONTAMINATE_ASSEMBLY as DECONTAMINATE_ASSEMBLY_CONTIG } from './workflows/decontaminate_assembly.nf'
include { DECONTAMINATE_ASSEMBLY as DECONTAMINATE_ASSEMBLY_SCAFFOLD } from './workflows/decontaminate_assembly.nf'
include { GENERATE_DECONTAM_EVIDENCE } from './workflows/generate_decontam_evidence.nf'


/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { BAM_TO_FASTQ } from './modules/bam_to_fastq.nf'
include { BUILD_MERYL_DB } from './modules/build_meryl_db.nf'
include { TRIM_HIC } from './modules/trim_hic.nf'
include { HIFIASM } from './modules/hifiasm.nf'

include { FIND_MITO_REFERENCE } from './modules/find_mito_reference.nf'
include { MITOHIFI } from './modules/mitohifi.nf'
include { FILTER_MITO_CONTIGS as FILTER_MITO_CONTIGS_HAP1 } from './modules/filter_mito_contigs.nf'
include { FILTER_MITO_CONTIGS as FILTER_MITO_CONTIGS_HAP2 } from './modules/filter_mito_contigs.nf'
include { MITO_CIRCULAR_MAP } from './modules/mito_circular_map.nf'

include { PURGE_DUPS } from './modules/purge_dups.nf'
include { CORRECT_MISASSEMBLIES as CORRECT_MISASSEMBLIES_CONTIG } from './modules/correct_misassemblies.nf'
include { CORRECT_MISASSEMBLIES as CORRECT_MISASSEMBLIES_SCAFFOLD } from './modules/correct_misassemblies.nf'
include { MAP_HIC_TO_ASSEMBLY } from './modules/map_hic_to_assembly.nf'
include { MAP_HIC_TO_ASSEMBLY as MAP_HIC_TO_SCAFFOLD } from './modules/map_hic_to_assembly.nf'
include { MAP_HIC_TO_ASSEMBLY as MAP_HIC_TO_FINAL } from './modules/map_hic_to_assembly.nf'
include { FILTER_HIC_BAM } from './modules/filter_hic_bam.nf'
include { FILTER_HIC_BAM as FILTER_HIC_BAM_SCAFFOLD } from './modules/filter_hic_bam.nf'
include { FILTER_HIC_BAM as FILTER_HIC_BAM_FINAL } from './modules/filter_hic_bam.nf'
include { SCAFFOLD_HIC as SCAFFOLD_HIC_ROUND1 } from './modules/scaffold_hic.nf'
include { SCAFFOLD_HIC as SCAFFOLD_HIC_ROUND2 } from './modules/scaffold_hic.nf'
include { GAP_FILLING } from './modules/gap_filling.nf'
include { TIDK_EXPLORE; TIDK_SEARCH; TIDK_PLOT; TIDK_SUMMARIZE; COLLECT_TIDK_RESULTS } from './modules/tidk.nf'
include { TELOCLIP_EXTEND; COLLECT_TELOCLIP_STATS } from './modules/teloclip.nf'
include { QUAST_FINAL } from './modules/quast.nf'
include { HIC_BAM_METRICS as HIC_BAM_METRICS_CONTIG; HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_CONTIG } from './modules/hic_mapping_metrics.nf'
include { HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_CONTIGSCAF } from './modules/hic_mapping_metrics.nf'
include { HIC_BAM_METRICS as HIC_BAM_METRICS_SCAFFOLD; HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_SCAFFOLD } from './modules/hic_mapping_metrics.nf'
include { HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_SCAFFOLDSCAF } from './modules/hic_mapping_metrics.nf'
include { HIC_BAM_METRICS as HIC_BAM_METRICS_FINAL; HIC_PAIRS_METRICS as HIC_PAIRS_METRICS_FINAL } from './modules/hic_mapping_metrics.nf'
include { COMPILE_FINAL_QC } from './modules/compile_final_qc.nf'
include { SNAIL_PLOT as SNAIL_PLOT_FINAL } from './modules/snail_plot.nf'
include { CONTACT_MAP as CONTACT_MAP_FINAL } from './modules/contact_map.nf'
include { SETUP_PAFR; PAIRWISE_ALIGNMENT; COLLECT_PAIRWISE_RESULTS } from './modules/pairwise_alignment.nf'
include { RIPARIAN_PLOT } from './modules/riparian_plot.nf'
//include { SCAN_TELOMERES; COLLECT_TELOMERE_RESULTS } from './modules/scan_telomeres.nf'
include { SUMMARY_REPORT } from './modules/summary_report'
include { DOWNLOAD_BUSCO_DB } from './modules/download_busco_db.nf'
include { COVERAGE_BOOK } from './modules/coverage_book.nf'
include { HIC_COMPARTMENTS } from './modules/hic_compartments.nf'
include { HIC_TADS } from './modules/hic_tads.nf'
include { ASSEMBLY_REPORT } from './modules/assemblyReport.nf'
include { FINALIZE_ASSEMBLY } from './modules/finalize_assembly.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Parse sample sheet and create input channel
    ch_input = parseSampleSheet(params.sample_sheet)

    /*
    ========================================================================================
        STEP 0: Setup Decontamination Databases (if enabled) & BUSCO Database
        Runs in parallel with BAM conversion and assembly
        Only executes if decontamination is requested
    ========================================================================================
    */
    if (params.decon.run_on_contigs || params.decon.run_on_scaffolds) {
        SETUP_DECONTAM_DBS()
        
        // Store outputs for later use
        ch_gxdb_dir = SETUP_DECONTAM_DBS.out.gxdb_dir
        ch_diamond_db = SETUP_DECONTAM_DBS.out.diamond_db
        ch_taxdump_dir = SETUP_DECONTAM_DBS.out.taxdump_dir
    }

    /*
    // Debug: Print all channel contents
    ch_input.view { sample_id, hifi_bam, hic_r1, hic_r2 ->
        """
        ========================================
        Sample ID  : ${sample_id}
        HiFi BAM   : ${hifi_bam}
        Hi-C R1    : ${hic_r1}
        Hi-C R2    : ${hic_r2}
        ========================================
        """
    }
    */
    
    DOWNLOAD_BUSCO_DB(params.busco_lineage)
    ch_busco_db = DOWNLOAD_BUSCO_DB.out.db

    /*
    ========================================================================================
        STEP 0b: Find Closest Mitochondrial Reference (NCBI)
        Runs at pipeline start — no dependencies on reads or assembly
    ========================================================================================
    */
    FIND_MITO_REFERENCE(params.mitohifi_species)

    /*
    ========================================================================================
        STEP 1: Convert BAM to FASTQ
    ========================================================================================
    */
    
    BAM_TO_FASTQ(
        ch_input.map { sample_id, hifi_bam, hic_r1, hic_r2 ->
            tuple(sample_id, hifi_bam)
        }
    )

    
    /*
    ========================================================================================
        STEP 2: Trim Hi-C Reads
    ========================================================================================
    */
    TRIM_HIC(
        ch_input.map { sample_id, hifi_bam, hic_r1, hic_r2 ->
            tuple(sample_id, hic_r1, hic_r2)
        }
    )

    /*
    ========================================================================================
        STEP 3: Combine HiFi FASTQ with trimmed Hi-C reads
    ========================================================================================
    */
    TRIM_HIC.out.trimmed_reads
        .join(BAM_TO_FASTQ.out)
        .map { sample_id, hic_r1_trim, hic_r2_trim, hifi_fastq ->
            tuple(sample_id, hifi_fastq, hic_r1_trim, hic_r2_trim)
        }
        .set { ch_fastq_all }

    /*
    ========================================================================================
        STEP 3b: Assemble Mitochondrial Genome (MitoHiFi)
        Runs concurrently with HIFIASM — both consume HiFi reads independently
        Depends on: BAM_TO_FASTQ + FIND_MITO_REFERENCE
    ========================================================================================
    */
    MITOHIFI(
        BAM_TO_FASTQ.out,
        FIND_MITO_REFERENCE.out.ref_fasta,
        FIND_MITO_REFERENCE.out.ref_gb
    )

    MITO_CIRCULAR_MAP(MITOHIFI.out.annotation)

    /*
    ========================================================================================
        STEP 4: Assemble with Hifiasm
    ========================================================================================
    */
    
    HIFIASM(ch_fastq_all)

    /*
    ========================================================================================
        STEP 4a: Remove Mitochondrial Contigs from Nuclear Assembly
        Runs after HIFIASM (+ optional purge_dups) and after MITOHIFI completes.
        Filters each haplotype independently using the assembled mitogenome.
        The filtered assemblies flow into all downstream steps (Inspector,
        decontamination, scaffolding, etc.)
    ========================================================================================
    */
    HIFIASM.out.assemblies
        .map { sample_id, hap1_fasta, hap2_fasta ->
            tuple(sample_id, hap1_fasta)
        }
        .combine(
            MITOHIFI.out.mitogenome.map { sid, mfa -> tuple(sid, mfa) },
            by: 0
        )
        .map { sample_id, hap1_fasta, mito_fa ->
            tuple("${sample_id}_hap1", hap1_fasta, mito_fa)
        }
        .set { ch_mito_filter_hap1_input }

    HIFIASM.out.assemblies
        .map { sample_id, hap1_fasta, hap2_fasta ->
            tuple(sample_id, hap2_fasta)
        }
        .combine(
            MITOHIFI.out.mitogenome.map { sid, mfa -> tuple(sid, mfa) },
            by: 0
        )
        .map { sample_id, hap2_fasta, mito_fa ->
            tuple("${sample_id}_hap2", hap2_fasta, mito_fa)
        }
        .set { ch_mito_filter_hap2_input }

    FILTER_MITO_CONTIGS_HAP1(ch_mito_filter_hap1_input)
    FILTER_MITO_CONTIGS_HAP2(ch_mito_filter_hap2_input)

    // Reconstruct sample-level tuple
    FILTER_MITO_CONTIGS_HAP1.out.filtered
        .mix(FILTER_MITO_CONTIGS_HAP2.out.filtered)
        .map { haplotype_id, fasta ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            def hap_num = haplotype_id.contains('_hap1') ? 1 : 2
            tuple(sample_id, hap_num, fasta)
        }
        .groupTuple(by: 0, size: 2)
        .map { sample_id, hap_nums, fastas ->
            def sorted = [hap_nums, fastas].transpose().sort { it[0] }
            tuple(sample_id, sorted[0][1], sorted[1][1])
        }
        .set { ch_mito_filtered }

    /*
    ========================================================================================
        STEP 4a: PURGE DUPLICATES (Optional)
        Remove haplotig duplications from HIFIASM output
        Only runs if params.run_purge_dups is true
        Runs BEFORE Inspector correction (if enabled)
    ========================================================================================
    */
    if (params.run_purge_dups) {
        ch_mito_filtered                              // <-- mito-filtered, not HIFIASM.out
            .flatMap { sample_id, hap1_fasta, hap2_fasta ->
                [
                    tuple("${sample_id}_hap1", hap1_fasta),
                    tuple("${sample_id}_hap2", hap2_fasta)
                ]
            }
            .map { haplotype_id, assembly ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(sample_id, haplotype_id, assembly)
            }
            .combine(BAM_TO_FASTQ.out, by: 0)
            .map { sample_id, haplotype_id, assembly, hifi_reads ->
                tuple(haplotype_id, assembly, hifi_reads)
            }
            .set { ch_purge_dups_input }

        PURGE_DUPS(ch_purge_dups_input)

        PURGE_DUPS.out.purged_assembly
            .map { haplotype_id, purged_fa ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = haplotype_id.contains('_hap1') ? 1 : 2
                tuple(sample_id, hap_num, purged_fa)
            }
            .groupTuple(by: 0, size: 2)
            .map { sample_id, hap_nums, fastas ->
                def sorted = [hap_nums, fastas].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_hifiasm_output }

    } else {
        ch_mito_filtered                              // <-- mito-filtered, not HIFIASM.out
            .set { ch_hifiasm_output }
    }


    /*
    ========================================================================================
        STEP 4.5: Optional Misassembly Correction of Contig Assemblies (Inspector)
        Runs before decontamination if enabled
        Uses HiFi reads to identify and break structural errors
        uses either purged or original assemblies based on params.run_purge_dups
    ========================================================================================
    */
    if (params.inspector_run_on_contigs) {
        // Split assemblies into individual haplotypes for correction
        ch_hifiasm_output
            .flatMap { sample_id, hap1_fasta, hap2_fasta ->
                [
                    tuple("${sample_id}_hap1", hap1_fasta),
                    tuple("${sample_id}_hap2", hap2_fasta)
                ]
            }
            .set { ch_contigs_for_correction }
        
        // Combine each haplotype with its HiFi reads for correction
        ch_contigs_for_correction
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(sample_id, haplotype_id, fasta)
            }
            .combine(BAM_TO_FASTQ.out, by: 0)
            .map { sample_id, haplotype_id, fasta, hifi_fastq ->
                // Create correction parameters map for CONTIGS
                def correction_params = [
                    min_depth: params.inspector_contig_min_depth,
                    min_contig_length: params.inspector_contig_min_contig_length,
                    min_contig_length_assemblyerror: params.inspector_contig_min_contig_length_assemblyerror,
                    min_assembly_error_size: params.inspector_contig_min_assembly_error_size,
                    max_assembly_error_size: params.inspector_contig_max_assembly_error_size,
                    skip_baseerror: params.inspector_contig_skip_baseerror
                ]
                tuple(haplotype_id, fasta, hifi_fastq, "contig", correction_params)
            }
            .set { ch_correction_input }
        
        // Run misassembly correction
        CORRECT_MISASSEMBLIES_CONTIG(ch_correction_input)
        
        // Use corrected assemblies for downstream processing
        ch_assemblies_for_decontam = CORRECT_MISASSEMBLIES_CONTIG.out.corrected
        
    } else {
        // Use original assemblies if correction not requested
        ch_hifiasm_output
            .flatMap { sample_id, hap1_fasta, hap2_fasta ->
                [
                    tuple("${sample_id}_hap1", hap1_fasta),
                    tuple("${sample_id}_hap2", hap2_fasta)
                ]
            }
            .set { ch_assemblies_for_decontam }
    }

    /*
    ========================================================================================
        STEP 5: Optional Decontamination of Contig Assemblies
        Runs in parallel with Hi-C mapping preparation
        Databases were already set up in STEP 0 (parallel with assembly)
        Works on either original contigs OR corrected contigs (if Inspector was run)
    ========================================================================================
    */
    if (params.decon.run_on_contigs) {
        // Run decontamination (parallel across all haplotypes)
        DECONTAMINATE_ASSEMBLY_CONTIG(
            ch_assemblies_for_decontam,
            ch_gxdb_dir,
            "contig"
        )
        
        
        // Use decontaminated assemblies for downstream Hi-C mapping
        // Need to add sample_id back for joining with Hi-C reads
        DECONTAMINATE_ASSEMBLY_CONTIG.out.decontaminated
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(haplotype_id, sample_id, fasta)
            }
            .set { ch_individual_haplotypes }
    } else {
        // Use corrected or original assemblies (depending on Inspector setting)
        ch_assemblies_for_decontam
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(haplotype_id, sample_id, fasta)
            }
            .set { ch_individual_haplotypes }
    }

    /*
    ========================================================================================
        STEP 6: Map Hi-C to Assemblies (contigs or decontaminated contigs)
    ========================================================================================
    */
    // Combine each haplotype with its corresponding trimmed Hi-C reads
    ch_individual_haplotypes
        .map { haplotype_id, sample_id, fasta ->
            tuple(sample_id, haplotype_id, fasta)
        }
        .combine(TRIM_HIC.out.trimmed_reads, by: 0)
        .map { sample_id, haplotype_id, fasta, hic_r1, hic_r2 ->
            tuple(haplotype_id, fasta, hic_r1, hic_r2, "contig")
        }
        .set { ch_hic_mapping_input }
    
    // Map Hi-C reads to assemblies
    MAP_HIC_TO_ASSEMBLY(ch_hic_mapping_input)

    // checkpoint: contig_raw_map (BAM-level only)
    MAP_HIC_TO_ASSEMBLY.out.bam
    .map { hap, stage, bam, bai -> tuple(hap, "contig_raw_map", bam, bai) }
    .set { ch_hic_raw_bam_for_qc }

    HIC_BAM_METRICS_CONTIG(ch_hic_raw_bam_for_qc)

    // Prepare assemblies channel for QC (haplotype_id, assembly_fasta)
    ch_individual_haplotypes
        .map { haplotype_id, sample_id, fasta ->
            tuple(haplotype_id, fasta)
        }
        .set { ch_assemblies_for_qc }


    /*
    ========================================================================================
        STEP 7: Filter Hi-C BAM Files
    ========================================================================================
    */
    
    // Combine BAM files with assemblies for filtering
    MAP_HIC_TO_ASSEMBLY.out.bam
        .join(ch_assemblies_for_qc)
        .map { haplotype_id, stage, bam, bai, assembly_fasta ->
            tuple(haplotype_id, stage, bam, bai, assembly_fasta)
        }
        .set { ch_bam_with_assembly }
    
    // Filter BAM files to remove invalid pairs and duplicates
    FILTER_HIC_BAM(ch_bam_with_assembly)

    // checkpoint: contig_filtered (pairs-level + retention)
    // After join: (hap, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats)
    // Note: No AGP at this checkpoint - pass empty list [] for optional AGP
    FILTER_HIC_BAM.out.pairs
        .join(FILTER_HIC_BAM.out.parse_stats)
        .join(FILTER_HIC_BAM.out.dedup_stats)
        .map { hap, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats ->
            // no agp at this checkpoint - pass empty list for optional path
            tuple(hap, "contig_filtered", pairs_gz, [], parse_stats, dedup_stats)
        }
        .set { ch_hic_pairs_contig_filtered_for_qc }

    HIC_PAIRS_METRICS_CONTIG(ch_hic_pairs_contig_filtered_for_qc)
    
    /*
    ========================================================================================
        STEP 8: First Round of Scaffolding with Hi-C
    ========================================================================================
    */
    
    // Prepare input for scaffolding: (haplotype_id, filtered_bam, bai, assembly_fasta, round, round_params)
    FILTER_HIC_BAM.out.bam
        .join(ch_assemblies_for_qc)
        .map { haplotype_id, stage, bam, bai, assembly_fasta ->
            tuple(haplotype_id, bam, bai, assembly_fasta, "round1", params.yahs_round1)
        }
        .set { ch_scaffolding_round1_input }
    
    // Run first round of Hi-C scaffolding
    SCAFFOLD_HIC_ROUND1(ch_scaffolding_round1_input)

    // checkpoint: scaffold_space (relabel contigs to scaffolds using AGP; no remapping)
    // After join: (hap, stage, pairs_gz, agp, stage2, parse_stats, stage3, dedup_stats)
    FILTER_HIC_BAM.out.pairs
        .join(SCAFFOLD_HIC_ROUND1.out.agp)
        .join(FILTER_HIC_BAM.out.parse_stats)
        .join(FILTER_HIC_BAM.out.dedup_stats)
        .map { hap, stage, pairs_gz, agp, stage2, parse_stats, stage3, dedup_stats ->
            tuple(hap, "scaffold_space", pairs_gz, agp, parse_stats, dedup_stats)
        }
        .set { ch_hic_pairs_scaffold_space_for_qc }

    HIC_PAIRS_METRICS_CONTIGSCAF(ch_hic_pairs_scaffold_space_for_qc)

    /*
    ========================================================================================
        STEP 8.5: Optional Misassembly Correction of Scaffolded Assemblies (Inspector)
        Runs after scaffolding, before scaffold decontamination if enabled
        Uses HiFi reads to identify and break structural errors in scaffolds
    ========================================================================================
    */
    if (params.inspector_run_on_scaffolds) {
        // Combine each scaffolded haplotype with its HiFi reads for correction
        SCAFFOLD_HIC_ROUND1.out.scaffolds
            .map { haplotype_id, scaffold ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(sample_id, haplotype_id, scaffold)
            }
            .combine(BAM_TO_FASTQ.out, by: 0)
            .map { sample_id, haplotype_id, scaffold, hifi_fastq ->
                // Create correction parameters map for SCAFFOLDS
                def correction_params = [
                    min_depth: params.inspector_scaffold_min_depth,
                    min_contig_length: params.inspector_scaffold_min_contig_length,
                    min_contig_length_assemblyerror: params.inspector_scaffold_min_contig_length_assemblyerror,
                    min_assembly_error_size: params.inspector_scaffold_min_assembly_error_size,
                    max_assembly_error_size: params.inspector_scaffold_max_assembly_error_size,
                    skip_baseerror: params.inspector_scaffold_skip_baseerror
                ]
                tuple(haplotype_id, scaffold, hifi_fastq, "scaffold", correction_params)
            }
            .set { ch_scaffold_correction_input }
        
        // Run misassembly correction on scaffolds
        CORRECT_MISASSEMBLIES_SCAFFOLD(ch_scaffold_correction_input)
        
        // Use corrected scaffolds for downstream processing
        ch_scaffolds_for_decontam = CORRECT_MISASSEMBLIES_SCAFFOLD.out.corrected
        
    } else {
        // Use original scaffolds if correction not requested
        ch_scaffolds_for_decontam = SCAFFOLD_HIC_ROUND1.out.scaffolds
    }

    /*
    ========================================================================================
        STEP 9: Optional Decontamination of Scaffolded Assemblies
        Runs after scaffolding (and optional correction) is complete
        Databases were already set up in STEP 0
        Works on either original scaffolds OR corrected scaffolds (if Inspector was run)
    ========================================================================================
    */
    if (params.decon.run_on_scaffolds) {
        // Decontaminate scaffolds (parallel across all haplotypes)
        DECONTAMINATE_ASSEMBLY_SCAFFOLD(
            ch_scaffolds_for_decontam,
            ch_gxdb_dir,
            "scaffold"
        )

        // Store final decontaminated scaffolds
        ch_final_scaffolds = DECONTAMINATE_ASSEMBLY_SCAFFOLD.out.decontaminated
    } else {
        // Use corrected or original scaffolds (depending on Inspector setting)
        ch_final_scaffolds = ch_scaffolds_for_decontam
    }
    
    /*
    ========================================================================================
        STEP 10-12: Conditional Second Round of Scaffolding
        Only runs if:
        - Inspector correction on scaffolds is enabled, OR
        - Decontamination on scaffolds is enabled, OR
        - Explicitly enabled via --run_scaffold_round2 true
        
        Can be explicitly disabled via --run_scaffold_round2 false
    ========================================================================================
    */
    if (params.run_scaffold_round2) {
        
        log.info "[INFO] Running second round of scaffolding (scaffold correction or decontamination was performed)"
        
        /*
        ========================================================================================
            STEP 10: Map Hi-C to Final Scaffolds (corrected/decontaminated if those options chosen)
        ========================================================================================
        */
        // Extract sample_id from haplotype_id and combine with trimmed Hi-C reads
        ch_final_scaffolds
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(sample_id, haplotype_id, fasta)
            }
            .combine(TRIM_HIC.out.trimmed_reads, by: 0)
            .map { sample_id, haplotype_id, fasta, hic_r1, hic_r2 ->
                tuple(haplotype_id, fasta, hic_r1, hic_r2, "scaffold")
            }
            .set { ch_hic_scaffold_mapping_input }
        
        // Map Hi-C reads to final scaffolds
        MAP_HIC_TO_SCAFFOLD(ch_hic_scaffold_mapping_input)

        // checkpoint: scaffold_round2_raw_map (BAM-level only)
        MAP_HIC_TO_SCAFFOLD.out.bam
            .map { hap, stage, bam, bai -> tuple(hap, "scaffold_round2_raw_map", bam, bai) }
            .set { ch_hic_scaffold_raw_bam_for_qc }

        HIC_BAM_METRICS_SCAFFOLD(ch_hic_scaffold_raw_bam_for_qc)

        /*
        ========================================================================================
            STEP 11: Filter Hi-C BAM Files mapped to final scaffolds
        ========================================================================================
        */
        
        // Combine BAM files with scaffolds for filtering
        MAP_HIC_TO_SCAFFOLD.out.bam
            .join(ch_final_scaffolds)
            .map { haplotype_id, stage, bam, bai, scaffold_fasta ->
                tuple(haplotype_id, stage, bam, bai, scaffold_fasta)
            }
            .set { ch_bam_with_scaffold }
        
        // Filter BAM files to remove invalid pairs and duplicates
        FILTER_HIC_BAM_SCAFFOLD(ch_bam_with_scaffold)

        // checkpoint: scaffold_round2_filtered (pairs-level + retention)
        // After join: (hap, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats)
        FILTER_HIC_BAM_SCAFFOLD.out.pairs
            .join(FILTER_HIC_BAM_SCAFFOLD.out.parse_stats)
            .join(FILTER_HIC_BAM_SCAFFOLD.out.dedup_stats)
            .map { hap, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats ->
                // already mapped to scaffold1 names, so no AGP needed here - pass empty list
                tuple(hap, "scaffold_round2_filtered", pairs_gz, [], parse_stats, dedup_stats)
            }
        .set { ch_hic_pairs_scaffold_round2_filtered_for_qc }

        HIC_PAIRS_METRICS_SCAFFOLD(ch_hic_pairs_scaffold_round2_filtered_for_qc)


        /*
        ========================================================================================
            STEP 12: Second Round of Scaffolding (Iterative Improvement)
            Uses filtered Hi-C BAM from corrected/decontaminated scaffolds
            This allows YaHS to potentially improve scaffolding based on the corrected structure
        ========================================================================================
        */
        
        // Prepare input for second scaffolding: (haplotype_id, filtered_bam, bai, scaffold_fasta, round, round_params)
        FILTER_HIC_BAM_SCAFFOLD.out.bam
            .join(ch_final_scaffolds)
            .map { haplotype_id, stage, bam, bai, scaffold_fasta ->
                tuple(haplotype_id, bam, bai, scaffold_fasta, "round2", params.yahs_round2)
            }
            .set { ch_second_scaffolding_input }
        
        // Run second round of Hi-C scaffolding (reusing same module with different round parameter)
        SCAFFOLD_HIC_ROUND2(ch_second_scaffolding_input)
        
        // Store final scaffolds from second round
        ch_final_scaffolds_round2 = SCAFFOLD_HIC_ROUND2.out.scaffolds

        // checkpoint: scaffold_round2_space (relabel scaffold1 -> scaffold2 using ROUND2 AGP; no remapping)
        // After join: (hap, stage, pairs_gz, agp, stage2, parse_stats, stage3, dedup_stats)
        FILTER_HIC_BAM_SCAFFOLD.out.pairs
            .join(SCAFFOLD_HIC_ROUND2.out.agp)
            .join(FILTER_HIC_BAM_SCAFFOLD.out.parse_stats)
            .join(FILTER_HIC_BAM_SCAFFOLD.out.dedup_stats)
            .map { hap, stage, pairs_gz, agp, stage2, parse_stats, stage3, dedup_stats ->
                tuple(hap, "scaffold_round2_space", pairs_gz, agp, parse_stats, dedup_stats)
            }
            .set { ch_hic_pairs_scaffold_round2_space_for_qc }

        HIC_PAIRS_METRICS_SCAFFOLDSCAF(ch_hic_pairs_scaffold_round2_space_for_qc)
        
    } else {
        log.info "[INFO] Skipping second round of scaffolding (no scaffold correction or decontamination)"
        
        // Create empty channel for round 2 scaffolds when not running
        ch_final_scaffolds_round2 = Channel.empty()
    }

    /*
    ========================================================================================
        STEP 13: Gap Filling
        Fills gaps in final scaffolded assemblies using HiFi reads
        Operates on the final scaffold output from either:
        - Round 2 scaffolding (if round 2 was run)
        - Decontaminated scaffolds (if decontamination on scaffolds was run)
        - Corrected scaffolds (if correction on scaffolds was run)
        - Original scaffolds (from round 1)
    ========================================================================================
    */
    
    // Determine which scaffolds to gap fill based on what was run
    if (params.run_scaffold_round2) {
        // Use round 2 scaffolds
        ch_scaffolds_for_gap_filling = ch_final_scaffolds_round2
    } else {
        // Use round 1 final scaffolds (corrected/decontaminated if those options were chosen)
        ch_scaffolds_for_gap_filling = ch_final_scaffolds
    }
    
    // Combine scaffolds with HiFi reads for gap filling
    ch_scaffolds_for_gap_filling
        .map { haplotype_id, scaffold ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            tuple(sample_id, haplotype_id, scaffold)
        }
        .combine(BAM_TO_FASTQ.out, by: 0)
        .map { sample_id, haplotype_id, scaffold, hifi_fastq ->
            tuple(haplotype_id, scaffold, hifi_fastq)
        }
        .set { ch_gap_filling_input }
    
    // Run gap filling
    GAP_FILLING(ch_gap_filling_input)

    /*
    ========================================================================================
        STEP 13b: Teloclip — Extend scaffolds with missing telomeres (Optional)
        Maps raw HiFi reads back to gap-filled scaffolds to find soft-clipped
        alignments at scaffold ends containing telomeric motifs, then appends
        the overhang sequence to recover missing telomeres.
    ========================================================================================
    */
    if (params.run_teloclip_extend) {
        // Combine gap-filled assemblies with HiFi reads (same join pattern as gap filling)
        GAP_FILLING.out.filled_assembly
            .map { haplotype_id, filled_fa ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(sample_id, haplotype_id, filled_fa)
            }
            .combine(BAM_TO_FASTQ.out, by: 0)
            .map { sample_id, haplotype_id, filled_fa, hifi_fastq ->
                tuple(haplotype_id, filled_fa, hifi_fastq)
            }
            .set { ch_teloclip_input }

        TELOCLIP_EXTEND(ch_teloclip_input)

        // Collect teloclip stats across all haplotypes
        COLLECT_TELOCLIP_STATS(
            TELOCLIP_EXTEND.out.stats.map { haplotype_id, stats -> stats }.collect()
        )

        // The teloclip-extended assembly becomes the "final" assembly
        ch_final_assembly = TELOCLIP_EXTEND.out.extended_assembly
        ch_teloclip_stats_for_report = COLLECT_TELOCLIP_STATS.out.stats

    } else {
        // No teloclip — gap-filled assembly IS the final assembly
        ch_final_assembly = GAP_FILLING.out.filled_assembly
        ch_teloclip_stats_for_report = Channel.of(file('NO_TELOCLIP'))
    }

    // =========================================================================
    //  FINALIZE ASSEMBLY — now uses ch_final_assembly (post-teloclip if enabled)
    // =========================================================================
    FINALIZE_ASSEMBLY(ch_final_assembly)

    /*
    ========================================================================================
        STEP 14: Finalization
    ========================================================================================
    */
    // 1. Contact Maps for Final Assemblies
    //    REPLACE: GAP_FILLING.out.filled_assembly → ch_final_assembly
    if (params.make_final_contact_maps) {
        ch_final_assembly
            .map { haplotype_id, final_fa ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                tuple(sample_id, haplotype_id, final_fa)
            }
            .combine(TRIM_HIC.out.trimmed_reads, by: 0)
            .map { sample_id, haplotype_id, final_fa, hic_r1, hic_r2 ->
                tuple(haplotype_id, final_fa, hic_r1, hic_r2, 'final')
            }
            .set { ch_final_hic_map_input }

        MAP_HIC_TO_FINAL(ch_final_hic_map_input)

        // (rest of contact map wiring is the same but uses ch_final_assembly
        //  instead of GAP_FILLING.out.filled_assembly for joins)

        MAP_HIC_TO_FINAL.out.bam
            .map { haplotype_id, stage, bam, bai -> tuple(haplotype_id, "final_raw_map", bam, bai) }
            .set { ch_final_raw_bam_qc }

        HIC_BAM_METRICS_FINAL(ch_final_raw_bam_qc)

        MAP_HIC_TO_FINAL.out.bam
            .join(ch_final_assembly, by: 0)
            .map { haplotype_id, stage, bam, bai, final_fa ->
                tuple(haplotype_id, "final", bam, bai, final_fa)
            }
            .set { ch_final_filter_input }

        FILTER_HIC_BAM_FINAL(ch_final_filter_input)

        FILTER_HIC_BAM_FINAL.out.pairs
            .join(FILTER_HIC_BAM_FINAL.out.parse_stats)
            .join(FILTER_HIC_BAM_FINAL.out.dedup_stats)
            .map { hap, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats ->
                tuple(hap, "final_filtered", pairs_gz, [], parse_stats, dedup_stats)
            }
            .set { ch_final_pairs_qc }

        HIC_PAIRS_METRICS_FINAL(ch_final_pairs_qc)

        FILTER_HIC_BAM_FINAL.out.pairs
            .join(ch_final_assembly)
            .map { haplotype_id, stage, pairs_gz, final_fasta ->
                tuple(haplotype_id, pairs_gz, final_fasta, "final")
            }
            .set { ch_contact_map_final_input }

        CONTACT_MAP_FINAL(ch_contact_map_final_input)

        if (params.hic_balance ) {
            HIC_COMPARTMENTS(
                CONTACT_MAP_FINAL.out.mcool,
                params.compartment_resolution ?: 250000,
                params.compartment_min_contig_bp ?: 5000000,
                params.compartment_max_contigs ?: 30
            )

            HIC_TADS(
                CONTACT_MAP_FINAL.out.mcool,
                params.tad_resolution ?: 50000,
                params.tad_window_bp ?: 500000,
                params.tad_min_contig_bp ?: 5000000,
                params.tad_max_contigs ?: 0
            )
        }
    }

    // 3. Dotplots of final assemblies vs each other (hap1 vs hap2)
    /*
    ========================================================================================
        PAIRWISE GENOME ALIGNMENTS
        Generates all pairwise alignments for dotplot visualization grid
    ========================================================================================
    */
    if (params.run_pairwise_alignments) {
        SETUP_PAFR()

        // Collect all final assemblies, SORT for deterministic ordering, then generate pairs
        ch_final_assembly
            .toSortedList { a, b -> a[0] <=> b[0] }
            .flatMap { assemblies ->
                def pairs = []

                if (params.pairwise_alignment_mode == 'within_sample') {
                    // Only compare hap1 vs hap2 within each sample
                    def grouped = assemblies.groupBy { haplotype_id, fasta ->
                        haplotype_id.replaceAll(/_hap[12]$/, '')
                    }
                    grouped.each { sample_id, haps ->
                        if (haps.size() == 2) {
                            def sorted = haps.sort { it[0] }
                            pairs << tuple(sorted[0][0], sorted[0][1], sorted[1][0], sorted[1][1])
                        }
                    }
                } else {
                    // Generate all unique pairs (i < j to avoid duplicates and self-alignments)
                    for (int i = 0; i < assemblies.size(); i++) {
                        for (int j = i + 1; j < assemblies.size(); j++) {
                            def (id1, fa1) = assemblies[i]
                            def (id2, fa2) = assemblies[j]
                            pairs << tuple(id1, fa1, id2, fa2)
                        }
                    }
                }

                return pairs
            }
            .set { ch_pairwise_input }

        // FIX: pass SETUP_PAFR.out.ready as second argument
        PAIRWISE_ALIGNMENT(ch_pairwise_input, SETUP_PAFR.out.ready)

        // Riparian plot input — uses ch_final_assembly
        ch_paf_with_asm1 = PAIRWISE_ALIGNMENT.out.paf
            .map { hap1, hap2, paf -> tuple(hap1, hap2, paf) }
            .combine(ch_final_assembly, by: 0)

        ch_riparian_input = ch_paf_with_asm1
            .map { hap1, hap2, paf, fasta1 -> tuple(hap2, hap1, fasta1, paf) }
            .combine(ch_final_assembly, by: 0)
            .map { hap2, hap1, fasta1, paf, fasta2 ->
                tuple(hap1, fasta1, hap2, fasta2, paf)
            }

        RIPARIAN_PLOT(ch_riparian_input)

        COLLECT_PAIRWISE_RESULTS(
            PAIRWISE_ALIGNMENT.out.qc.map { id1, id2, qc_file -> qc_file }.collect()
        )
        ch_pairwise_summary = COLLECT_PAIRWISE_RESULTS.out.summary
    } else {
        ch_pairwise_summary = Channel.of(file('NO_PAIRWISE'))
    }

    /*
    ========================================================================================
        QUAST Final Comparison - All Gap-Filled Genomes
    ========================================================================================
    */
    // Collect all gap-filled assemblies with their labels for cross-sample comparison
    ch_final_assembly
        .map { haplotype_id, assembly -> tuple(haplotype_id, assembly) }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { list ->
            def labels = list.collect { it[0] }
            def assemblies = list.collect { it[1] }
            tuple(assemblies, labels)
        }
        .set { ch_quast_final_input }

    ch_quast_final_input
        .map { assemblies, labels -> assemblies }
        .set { ch_quast_assemblies }

    ch_quast_final_input
        .map { assemblies, labels -> labels }
        .set { ch_quast_labels }

    QUAST_FINAL(
        ch_quast_assemblies.collect(),
        ch_quast_labels.flatten().collect()
    )

    /*
    ========================================================================================
        Telomere Detection in Final Assemblies
    ========================================================================================
    */
    // 1. Explore: discover telomere motif de novo
    TIDK_EXPLORE(ch_final_assembly)

    // 2. Search: windowed telomere repeat quantification
    TIDK_SEARCH(ch_final_assembly)

    // 3. Plot: SVG visualization per haplotype
    TIDK_PLOT(TIDK_SEARCH.out.search_tsv)

    // 4. Summarize: per-scaffold presence/absence (backward-compatible with old format)
    TIDK_SUMMARIZE(TIDK_SEARCH.out.search_tsv)

    // 5. Collect all results
    COLLECT_TIDK_RESULTS(
        TIDK_SUMMARIZE.out.summary.map   { haplotype_id, summary   -> summary   }.collect(),
        TIDK_SUMMARIZE.out.telomeres.map { haplotype_id, telomeres -> telomeres }.collect()
    )

    // 6. NCBI output files for GenBank submission (if enabled)
    
    /*
    ========================================================================================
        QC Steps
    ========================================================================================
    */

    /*
    ========================================================================================
        BUILD MERYL DATABASE (Once per sample)
        Runs in parallel with trimming and assembly
        Reused across ALL assembly QC steps for dramatic speedup
    ========================================================================================
    */
    BUILD_MERYL_DB(BAM_TO_FASTQ.out)

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
    HIC_QC_RAW(
        ch_input.map { sample_id, hifi_bam, hic_r1, hic_r2 ->
            tuple(sample_id, hic_r1, hic_r2)
        },
        "raw"
    )

    /*
    ========================================================================================
        QC HiFi Reads
    ========================================================================================
    */
    HIFI_QC(
        BAM_TO_FASTQ.out
    )

    /*
    ========================================================================================
        QC Trimmed Hi-C Reads
    ========================================================================================
    */
    HIC_QC_TRIMMED(
        TRIM_HIC.out.trimmed_reads,
        "trimmed"
    )
    
    /*
    ========================================================================================
        Assembly QC
    ========================================================================================
    */

    /*
    ========================================================================================
        QC Contig Genomes Assemblies
    ========================================================================================
    */
    ASSEMBLY_QC_INITIAL(
        HIFIASM.out.assemblies,
        BAM_TO_FASTQ.out,
        BUILD_MERYL_DB.out.meryl_db,
        ch_busco_db,
        'contig'
    )

    /*
    ========================================================================================
        QC Mito-Filtered Contig Assemblies
    ========================================================================================
    */
    ASSEMBLY_QC_MITO_FILTERED(
        ch_mito_filtered,
        BAM_TO_FASTQ.out,
        BUILD_MERYL_DB.out.meryl_db,
        ch_busco_db,
        'contig_mito_filtered'
    )

    /*
    ========================================================================================
        QC Purged Contig Genomes Assemblies
    ========================================================================================
    */
    if (params.run_purge_dups) {
        ASSEMBLY_QC_PURGED(
            ch_hifiasm_output,    // tuple(sample_id, hap1_fasta, hap2_fasta)
            BAM_TO_FASTQ.out,     // tuple(sample_id, hifi_fastq)
            BUILD_MERYL_DB.out.meryl_db,        // tuple(sample_id, meryl_db)
            ch_busco_db,          // path
            "contig_purged"       // string
        )
    }
    /*
    ========================================================================================
        QC Corrected Contig Genomes Assemblies (if Inspector was run)
    ========================================================================================
    */
    if (params.inspector_run_on_contigs) {
        // Re-pair corrected haplotypes by sample for QC
        CORRECT_MISASSEMBLIES_CONTIG.out.corrected
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
                tuple(sample_id, hap_num, fasta)
            }
            .groupTuple(by: 0, size: 2)
            .map { sample_id, hap_nums, fastas ->
                // Sort by haplotype number to ensure hap1, hap2 order
                def sorted = [hap_nums, fastas].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_corrected_paired }

        ASSEMBLY_QC_CONTIG_CORRECTED(
            ch_corrected_paired,
            BAM_TO_FASTQ.out,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'contig_corrected'
        )
    }

    /*
    ========================================================================================
        QC Contig Decontaminated Genome Assemblies
    ========================================================================================
    */
    if (params.decon.run_on_contigs) {
        // QC on decontaminated assemblies
        // Re-pair haplotypes by sample for QC
        DECONTAMINATE_ASSEMBLY_CONTIG.out.decontaminated
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
                tuple(sample_id, hap_num, fasta)
            }
            .groupTuple(by: 0, size: 2)
            .map { sample_id, hap_nums, fastas ->
                // Sort by haplotype number to ensure hap1, hap2 order
                def sorted = [hap_nums, fastas].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_decontam_paired }

        ASSEMBLY_QC_CONTIG_DECONTAM(
            ch_decontam_paired,
            BAM_TO_FASTQ.out,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'contig_decontam'
        )
    }

    /*
    ========================================================================================
        QC Scaffolded Genomes (Round 1)
    ========================================================================================
    */
    // Re-pair scaffolded haplotypes by sample
    SCAFFOLD_HIC_ROUND1.out.scaffolds
        .map { haplotype_id, scaffold ->
            // Extract sample_id and haplotype number
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
            tuple(sample_id, hap_num, scaffold)
        }
        .groupTuple(by: 0, size: 2)
        .map { sample_id, hap_nums, scaffolds ->
            // Sort by haplotype number to ensure hap1, hap2 order
            def sorted = [hap_nums, scaffolds].transpose().sort { it[0] }
            tuple(sample_id, sorted[0][1], sorted[1][1])
        }
        .set { ch_scaffolds_paired }

    ASSEMBLY_QC_SCAFFOLD(
        ch_scaffolds_paired,
        BAM_TO_FASTQ.out,
        BUILD_MERYL_DB.out.meryl_db,
        ch_busco_db,
        'scaffold'
    )

    /*
    ========================================================================================
        QC Corrected Scaffold Genomes Assemblies (if Inspector was run on scaffolds)
    ========================================================================================
    */
    if (params.inspector_run_on_scaffolds) {
        // Re-pair corrected scaffolds by sample for QC
        CORRECT_MISASSEMBLIES_SCAFFOLD.out.corrected
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
                tuple(sample_id, hap_num, fasta)
            }
            .groupTuple(by: 0, size: 2)
            .map { sample_id, hap_nums, fastas ->
                // Sort by haplotype number to ensure hap1, hap2 order
                def sorted = [hap_nums, fastas].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_scaffold_corrected_paired }

        ASSEMBLY_QC_SCAFFOLD_CORRECTED(
            ch_scaffold_corrected_paired,
            BAM_TO_FASTQ.out,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'scaffold_corrected'
        )
    }

    /*
    ========================================================================================
        QC Scaffold Decontaminated Genome Assemblies
    ========================================================================================
    */
    if (params.decon.run_on_scaffolds) {
        // QC on decontaminated assemblies
        // Re-pair haplotypes by sample for QC
        DECONTAMINATE_ASSEMBLY_SCAFFOLD.out.decontaminated
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
                tuple(sample_id, hap_num, fasta)
            }
            .groupTuple(by: 0, size: 2)
            .map { sample_id, hap_nums, fastas ->
                // Sort by haplotype number to ensure hap1, hap2 order
                def sorted = [hap_nums, fastas].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_decontam_scaffold_paired }

        ASSEMBLY_QC_SCAFFOLD_DECONTAM(
            ch_decontam_scaffold_paired,
            BAM_TO_FASTQ.out,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'scaffold_decontam'
        )
    }

    /*
    ========================================================================================
        QC Second Round Scaffolded Genomes (Conditional - only if round 2 was run)
    ========================================================================================
    */
    if (params.run_scaffold_round2) {
        // Re-pair second-round scaffolded haplotypes by sample
        ch_final_scaffolds_round2
            .map { haplotype_id, scaffold ->
                // Extract sample_id and haplotype number
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
                tuple(sample_id, hap_num, scaffold)
            }
            .groupTuple(by: 0, size: 2)
            .map { sample_id, hap_nums, scaffolds ->
                // Sort by haplotype number to ensure hap1, hap2 order
                def sorted = [hap_nums, scaffolds].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_scaffolds_round2_paired }

        ASSEMBLY_QC_SCAFFOLD_ROUND2(
            ch_scaffolds_round2_paired,
            BAM_TO_FASTQ.out,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'scaffold_round2'
        )
    }

    /*
    ========================================================================================
        QC Gap-Filled Genomes
    ========================================================================================
    */
    // Re-pair gap-filled haplotypes by sample
    GAP_FILLING.out.filled_assembly
        .map { haplotype_id, assembly ->
            // Extract sample_id and haplotype number
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
            tuple(sample_id, hap_num, assembly)
        }
        .groupTuple(by: 0, size: 2)
        .map { sample_id, hap_nums, assemblies ->
            // Sort by haplotype number to ensure hap1, hap2 order
            def sorted = [hap_nums, assemblies].transpose().sort { it[0] }
            tuple(sample_id, sorted[0][1], sorted[1][1])
        }
        .set { ch_gap_filled_paired }

    ASSEMBLY_QC_GAP_FILLED(
        ch_gap_filled_paired,
        BAM_TO_FASTQ.out,
        BUILD_MERYL_DB.out.meryl_db,
        ch_busco_db,
        'gap_filled'
    )

    /*
    ========================================================================================
        Generate Optional Decontamination Evidence
    ========================================================================================
    */
    if (params.decon.run_on_contigs) {
        // Optional: Generate evidence for decontamination decisions
        // This runs in parallel with Hi-C mapping preparation
        if (params.decon.make_blobtools_evidence) {
            GENERATE_DECONTAM_EVIDENCE(
                DECONTAMINATE_ASSEMBLY_CONTIG.out.decontaminated,
                DECONTAMINATE_ASSEMBLY_CONTIG.out.contaminants,
                DECONTAMINATE_ASSEMBLY_CONTIG.out.action_report,
                DECONTAMINATE_ASSEMBLY_CONTIG.out.taxonomy_report,
                BAM_TO_FASTQ.out,
                ch_diamond_db,
                ch_taxdump_dir
            )
        }
    }

    if (params.decon.run_on_scaffolds) {
        // Optional: Generate evidence for scaffold decontamination
        // This runs in parallel with scaffold QC
        if (params.decon.make_blobtools_evidence) {
            GENERATE_DECONTAM_EVIDENCE(
                DECONTAMINATE_ASSEMBLY_SCAFFOLD.out.decontaminated,
                DECONTAMINATE_ASSEMBLY_SCAFFOLD.out.contaminants,
                DECONTAMINATE_ASSEMBLY_SCAFFOLD.out.action_report,
                DECONTAMINATE_ASSEMBLY_SCAFFOLD.out.taxonomy_report,
                BAM_TO_FASTQ.out,
                ch_diamond_db,
                ch_taxdump_dir
            )
        }
    }

    /*
    ========================================================================================
        FINAL QC COMPILATION
        Collects all assembly QC summaries and Hi-C metrics into a single report
    ========================================================================================
    */
    
    // Collect all assembly QC summaries
    // Start with the ones that always run
    ch_all_assembly_summaries = ASSEMBLY_QC_INITIAL.out.assembly_summary
        .mix(ASSEMBLY_QC_MITO_FILTERED.out.assembly_summary)
        .mix(ASSEMBLY_QC_SCAFFOLD.out.assembly_summary)
        .mix(ASSEMBLY_QC_GAP_FILLED.out.assembly_summary)
    
    // Add conditional assembly QC outputs
    if (params.run_purge_dups) {
        ch_all_assembly_summaries = ch_all_assembly_summaries
            .mix(ASSEMBLY_QC_PURGED.out.assembly_summary)
    }

    if (params.inspector_run_on_contigs) {
        ch_all_assembly_summaries = ch_all_assembly_summaries
            .mix(ASSEMBLY_QC_CONTIG_CORRECTED.out.assembly_summary)
    }
    
    if (params.decon.run_on_contigs) {
        ch_all_assembly_summaries = ch_all_assembly_summaries
            .mix(ASSEMBLY_QC_CONTIG_DECONTAM.out.assembly_summary)
    }
    
    if (params.inspector_run_on_scaffolds) {
        ch_all_assembly_summaries = ch_all_assembly_summaries
            .mix(ASSEMBLY_QC_SCAFFOLD_CORRECTED.out.assembly_summary)
    }
    
    if (params.decon.run_on_scaffolds) {
        ch_all_assembly_summaries = ch_all_assembly_summaries
            .mix(ASSEMBLY_QC_SCAFFOLD_DECONTAM.out.assembly_summary)
    }
    
    if (params.run_scaffold_round2) {
        ch_all_assembly_summaries = ch_all_assembly_summaries
            .mix(ASSEMBLY_QC_SCAFFOLD_ROUND2.out.assembly_summary)
    }
    
    if (params.run_teloclip_extend) {
        ch_final_assembly
            .map { haplotype_id, fasta ->
                def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
                def hap_num = (haplotype_id =~ /_hap([12])$/)[0][1]
                tuple(sample_id, hap_num, fasta)
            }
            .groupTuple(by: 0, size: 2)
            .map { sample_id, hap_nums, fastas ->
                def sorted = [hap_nums, fastas].transpose().sort { it[0] }
                tuple(sample_id, sorted[0][1], sorted[1][1])
            }
            .set { ch_teloclip_assemblies_paired }

        ASSEMBLY_QC_TELOCLIP(
            ch_teloclip_assemblies_paired,
            BAM_TO_FASTQ.out,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'teloclip_extended'
        )

        ch_all_assembly_summaries = ch_all_assembly_summaries
            .mix(ASSEMBLY_QC_TELOCLIP.out.assembly_summary)

        ch_final_busco = ASSEMBLY_QC_TELOCLIP.out.busco_results
    } else {
        ch_final_busco = ASSEMBLY_QC_GAP_FILLED.out.busco_results
    }

    // Collect all BAM metrics
    // Start with ones that always run
    ch_all_bam_metrics = HIC_BAM_METRICS_CONTIG.out.metrics
    
    // Add conditional BAM metrics
    if (params.run_scaffold_round2) {
        ch_all_bam_metrics = ch_all_bam_metrics
            .mix(HIC_BAM_METRICS_SCAFFOLD.out.metrics)
    }
    
    if (params.make_final_contact_maps) {
        ch_all_bam_metrics = ch_all_bam_metrics
            .mix(HIC_BAM_METRICS_FINAL.out.metrics)
    }
    
    // Collect all pairs metrics
    // Start with ones that always run
    ch_all_pairs_metrics = HIC_PAIRS_METRICS_CONTIG.out.metrics
        .mix(HIC_PAIRS_METRICS_CONTIGSCAF.out.metrics)
    
    // Add conditional pairs metrics
    if (params.run_scaffold_round2) {
        ch_all_pairs_metrics = ch_all_pairs_metrics
            .mix(HIC_PAIRS_METRICS_SCAFFOLD.out.metrics)
            .mix(HIC_PAIRS_METRICS_SCAFFOLDSCAF.out.metrics)
    }
    
    if (params.make_final_contact_maps) {
        ch_all_pairs_metrics = ch_all_pairs_metrics
            .mix(HIC_PAIRS_METRICS_FINAL.out.metrics)
    }
    
    // Compile final QC report
    // Extract just the TSV files from tuples and collect
    COMPILE_FINAL_QC(
        ch_all_assembly_summaries.map { sample_id, qc_label, tsv -> tsv }.collect(),
        ch_all_bam_metrics.map { haplotype_id, checkpoint, tsv -> tsv }.collect(),
        ch_all_pairs_metrics.map { haplotype_id, checkpoint, tsv -> tsv }.collect()
    )

    // Make interactive HTML assembly viewer
    ASSEMBLY_REPORT(COMPILE_FINAL_QC.out.metrics)

    /*
    ========================================================================================
    SNAIL PLOTS FOR FINAL ASSEMBLIES
    ========================================================================================
    */
    // Join gap-filled assemblies with their BUSCO results
    // BUSCO output from ASSEMBLY_QC_TELOCLIP is per-haplotype
    ch_final_assembly
        .join(ch_final_busco)
        .map { haplotype_id, assembly, busco_dir ->
            tuple(haplotype_id, assembly, busco_dir, "final")
        }
        .set { ch_snail_plot_final_input }

    SNAIL_PLOT_FINAL(ch_snail_plot_final_input)
    /*
    ========================================================================================
        COVERAGE BOOK - HiFi coverage visualization for final assemblies
    ========================================================================================
    */

    // Combine gap-filled assemblies with HiFi reads for coverage book
    ch_final_assembly
        .map { haplotype_id, final_fa ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            tuple(sample_id, haplotype_id, final_fa)
        }
        .combine(BAM_TO_FASTQ.out, by: 0)
        .map { sample_id, haplotype_id, final_fa, hifi_fastq ->
            def meta = [
                id: haplotype_id,
                sample: sample_id
            ]
            tuple(meta, final_fa, hifi_fastq)
        }
        .set { ch_coverage_book_input }

    COVERAGE_BOOK(ch_coverage_book_input)

    // =========================================================================
    //  SUMMARY REPORT — Build manifest and call the process
    //
    //  Verified publishDir paths from each module:
    //    FINALIZE_ASSEMBLY  → ${params.outdir}/assembly/final
    //    GAP_FILLING        → ${params.outdir}/assembly/scaffold/gap_filling
    //    SNAIL_PLOT         → ${params.outdir}/snail_plots
    //    CONTACT_MAP        → ${params.outdir}/contact_maps
    //    PAIRWISE_ALIGNMENT → ${params.outdir}/pairwise_alignments
    //    COMPILE_FINAL_QC   → ${params.outdir}/qc/assembly
    //    ASSEMBLY_REPORT    → ${params.outdir}/reports
    // =========================================================================

    // ---- Final genome assemblies ----
    // FINALIZE_ASSEMBLY.out.assembly: tuple(haplotype_id, fasta)
    // publishDir: ${params.outdir}/assembly/final
    ch_manifest_assemblies = FINALIZE_ASSEMBLY.out.assembly
        .map { hap_id, fasta ->
            "assembly\t${hap_id}\t.\t${fasta.name}\tassembly/final"
        }

    // ---- Snail plots (final) ----
    // SNAIL_PLOT_FINAL.out.snail: tuple(haplotype_id, qc_label, svg)
    // publishDir: ${params.outdir}/snail_plots
    ch_manifest_snails = SNAIL_PLOT_FINAL.out.snail
        .map { hap_id, qc_label, svg ->
            "snail\t${hap_id}\t.\t${svg.name}\tsnail_plots"
        }

    // ---- Contact maps (conditional) ----
    // CONTACT_MAP_FINAL.out.contact_maps: tuple(haplotype_id, stage, png_files)
    // publishDir: ${params.outdir}/contact_maps
    if (params.make_final_contact_maps) {
        ch_manifest_contact_maps = CONTACT_MAP_FINAL.out.contact_maps
            .flatMap { hap_id, stage, pngs ->
                def png_list = pngs instanceof List ? pngs : [pngs]
                png_list.collect { png ->
                    "contact_map\t${hap_id}\t.\t${png.name}\tcontact_maps"
                }
            }
    } else {
        ch_manifest_contact_maps = Channel.empty()
    }

    // ---- Pairwise dotplots (conditional) ----
    // PAIRWISE_ALIGNMENT.out.dotplot: tuple(id1, id2, png)
    // publishDir: ${params.outdir}/pairwise_alignments
    if (params.run_pairwise_alignments) {
        ch_manifest_dotplots = PAIRWISE_ALIGNMENT.out.dotplot
            .map { id1, id2, png ->
                "dotplot\t${id1}\t${id2}\t${png.name}\tpairwise_alignments"
            }
    } else {
        ch_manifest_dotplots = Channel.empty()
    }

    // ---- Pairwise riparian plots (conditional) ----
    // RIPARIAN_PLOT.out.riparian: tuple(id1, id2, png)
    // publishDir: ${params.outdir}/pairwise_alignments
    if (params.run_pairwise_alignments) {
        ch_manifest_riparian = RIPARIAN_PLOT.out.riparian
            .map { id1, id2, png ->
                "riparian\t${id1}\t${id2}\t${png.name}\tpairwise_alignments"
            }
    } else {
        ch_manifest_riparian = Channel.empty()
    }

    // ---- tidk plot SVGs for manifest ----
    ch_manifest_tidk_plots = TIDK_PLOT.out.plot
        .map { hap_id, svg ->
            "tidk_plot\t${hap_id}\t.\t${svg.name}\ttelomeres/plots"
        }
    

    // ---- Compiled QC CSV from COMPILE_FINAL_QC ----
    // publishDir: ${params.outdir}/qc/assembly
    ch_manifest_compiled_csv = COMPILE_FINAL_QC.out.metrics
        .map { csv ->
            "compiled_qc\t.\t.\t${csv.name}\tqc/assembly"
        }

    // ---- QC trend plots (PNGs) from COMPILE_FINAL_QC ----
    // publishDir: ${params.outdir}/qc/assembly
    ch_manifest_qc_plots = COMPILE_FINAL_QC.out.plots
        .flatten()
        .map { png ->
            "qc_plot\t.\t.\t${png.name}\tqc/assembly"
        }

    // ---- Interactive HTML report from ASSEMBLY_REPORT ----
    // publishDir: ${params.outdir}/reports
    ch_manifest_assembly_report = ASSEMBLY_REPORT.out.report_html
        .map { html ->
            "assembly_report_html\t.\t.\t${html.name}\treports"
        }

    // ---- Mitogenome assembly ----
    // publishDir: ${params.outdir}/mitogenome/${sample_id}
    ch_manifest_mitogenome = MITOHIFI.out.mitogenome
        .map { sample_id, fasta ->
            "mitogenome\t${sample_id}\t.\t${fasta.name}\tmitogenome"
        }

    ch_manifest_mito_stats = MITOHIFI.out.stats
        .map { sample_id, tsv ->
            "mito_stats\t${sample_id}\t.\t${tsv.name}\tmitogenome"
        }

    ch_manifest_mito_circular = MITO_CIRCULAR_MAP.out.circular_map
        .map { sample_id, png ->
            "mito_gene_map\t${sample_id}\t.\t${png.name}\t/mitogenome"
        }

    // ---- Combine all manifest entries into a single TSV ----
    ch_manifest_assemblies
        .mix(ch_manifest_snails)
        .mix(ch_manifest_contact_maps)
        .mix(ch_manifest_dotplots)
        .mix(ch_manifest_riparian)
        .mix(ch_manifest_compiled_csv)
        .mix(ch_manifest_qc_plots)
        .mix(ch_manifest_assembly_report)
        .mix(ch_manifest_mitogenome)
        .mix(ch_manifest_mito_stats)
        .mix(ch_manifest_mito_circular)
        .mix(ch_manifest_tidk_plots)
        .collectFile(
            name: 'report_manifest.tsv',
            seed: 'type\tid\tid2\tfilename\tsubdir',
            newLine: true
        )
        .set { ch_report_manifest }

    // ---- Handle optional telomere output ----
    ch_telomere_for_report = COLLECT_TIDK_RESULTS.out.summary
        .ifEmpty(file('NO_TELOMERES'))

    // Collect mitogenome stats for report
    ch_mito_stats_for_report = MITOHIFI.out.stats
        .map { sample_id, tsv -> tsv }
        .collectFile(
            name: 'all_mito_stats.tsv',
            keepHeader: true,
            skip: 1
        )
        .ifEmpty(file('NO_MITO_STATS'))

    // ---- Call SUMMARY_REPORT ----
    SUMMARY_REPORT(
        ch_report_manifest,
        COMPILE_FINAL_QC.out.metrics,
        ch_telomere_for_report,
        ch_pairwise_summary,
        ch_mito_stats_for_report,
        ch_teloclip_stats_for_report       
    )
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