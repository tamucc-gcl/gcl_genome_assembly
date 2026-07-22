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
params.qc_mode = 'final_only'   // 'all_stages' (QC every checkpoint) | 'final_only'
params.species = null   // optional global organism-name override; else derived from taxid

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
params.mitohifi_genetic_code   = 2           // NCBI genetic code: 2=vertebrate mito, 5=invertebrate mito
params.mitohifi_ref_min_length = 14000       // Min reference mitogenome length for findMitoReference.py
params.mitohifi_perc_identity  = 50          // Min percent identity for read filtering
params.mitohifi_cov_cutoff     = 'auto'      // Coverage cutoff (auto or integer)
params.mitohifi_bloom_filter   = false       // Use bloom filter for large read sets
params.mitohifi_circular_map_dpi = 300

// Mito contig filtering (applied to HIFIASM output)
params.mitohifi_filter_min_identity = 90     // Min alignment identity to call a contig mitochondrial
params.mitohifi_filter_min_coverage = 50     // Min query coverage to call a contig mitochondrial

// ---- SPAdes (short-read assembly) ----
params.spades_kmers        = '21,33,55,77'   // -k
params.spades_cov_cutoff   = 'auto'          // --cov-cutoff (a number, or 'auto')
params.spades_mode         = '--isolate'     // '--isolate' | '--careful' | '--sc' | '' for none
params.spades_extra        = ''              // any extra spades.py flags, verbatim
params.spades_output_level = 'scaffolds'       // 'contigs' | 'scaffolds' — which FASTA feeds downstream

// ---- Genome-size estimation (jellyfish + GenomeScope2) ----
params.genomescope_kmer    = 21              // jellyfish -m / genomescope -k
params.jellyfish_hash_size = '5G'            // jellyfish count -s

// Purge duplicates parameters
params.run_purge_dups = false  // Default off - enable with --run_purge_dups true

// Redundans parameters
params.redundans_run_reduction      = true
params.redundans_run_scaffolding    = true
params.redundans_run_gapclosing     = true
params.redundans_identity           = 0.51
params.redundans_overlap            = 0.80
params.redundans_min_length         = 200
params.redundans_joins              = 5
params.redundans_linkratio          = 0.7
params.redundans_limit              = 0.2
params.redundans_mapq               = 10
params.redundans_iters              = 2
params.redundans_preset             = 'asm5'
params.redundans_index              = '4G'
params.redundans_minimap2reduce     = false
params.redundans_minimap2scaffold   = false
params.redundans_usebwa             = false
params.redundans_populate_scaffolds = false
params.redundans_norearrangements   = false
params.redundans_extra              = ''

// Shortread trim Parameters
params.shortread_trim                 = true      // raw + fastp; false = pass pre-cleaned reads through
params.shortread_fastp_cut_tail_quality = 20
params.shortread_fastp_cut_tail_window  = 4
params.shortread_fastp_length_required  = 30
params.shortread_fastp_extra          = ''

// Pilon parameters
params.run_pilon                      = false     // optional short-read polish
params.pilon_fix                      = 'all'
params.pilon_rounds                   = 1
params.pilon_extra                    = ''

// Assembly QC parameters
params.busco_lineage = 'eukaryota_odb10'   // fallback only; primary lineage is per-sample from taxonomy
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
params.finalize_min_scaffold_size = 1000000

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
params.taxid = null    // optional GLOBAL taxid for all samples; per-sample `taxid` column overrides; decontam falls back to params.decon_source_taxid
params.decon_source_taxid = null          // set your species taxid when possible

// Optional: adapter/vector removal
params.decon_run_fcs_adaptor = false      // Requires FCS-adaptor containers; enable only if configured
params.decon_fcsadaptor_mode = 'euk'      // 'euk' or 'prok'
params.decon_container_engine = 'singularity'

// Optional: generate evidence (coverage + taxonomy + blobtools plots)
params.decon_make_blobtools_evidence = false

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
include { CONTIG_ASSEMBLY } from './workflows/contig_assembly.nf'

//Organelle assembly/annotation
include { ORGANELLE_ASSEMBLY } from './workflows/organelle_assembly.nf'

//Short-read qc #Hi-C, SSL, TellSeq
include { SHORTREAD_QC as SHORTREAD_QC_RAW }     from './workflows/shortread_qc.nf'
include { SHORTREAD_QC as SHORTREAD_QC_TRIMMED } from './workflows/shortread_qc.nf'

// Assembly QC
include { ASSEMBLY_QC as ASSEMBLY_QC_INITIAL } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_MITO_FILTERED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_PURGED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_REDUNDANS } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_CONTIG_CORRECTED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_CONTIG_DECONTAM } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_CORRECTED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_DECONTAM } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_ROUND2 } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_GAP_FILLED } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_TELOCLIP } from './workflows/assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_FINAL } from './workflows/assembly_qc.nf'

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
include { RESOLVE_TAXONOMY } from './modules/resolve_taxonomy.nf'
include { buscoLineageFor; kingdomFlag; organismName } from './functions/taxonomy.nf'
include { DOWNLOAD_TAXDUMP } from './modules/download_taxdump.nf'

include { BAM_TO_FASTQ } from './modules/bam_to_fastq.nf'
include { BUILD_MERYL_DB } from './modules/build_meryl_db.nf'
include { TRIM_HIC } from './modules/trim_hic.nf'
include { ESTIMATE_GENOME_SIZE } from './modules/estimate_genome_size.nf'

include { FIND_MITO_REFERENCE } from './modules/find_mito_reference.nf'
include { FILTER_MITO_CONTIGS } from './modules/filter_mito_contigs.nf'

include { TRIM_SHORTREAD } from './modules/trim_shortread.nf'
include { REDUNDANS }      from './modules/redundans.nf'
include { PILON }          from './modules/pilon.nf'

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

// ── helper scripts declared as inputs so edits invalidate the cache ──
ch_compile_qc_script      = file("${projectDir}/r_scripts/compile_qc.R",                         checkIfExists: true)
ch_summary_report_script  = file("${projectDir}/r_scripts/generate_summary_report.R",            checkIfExists: true)
ch_dotplot_script         = file("${projectDir}/r_scripts/dotplot_paf.R",                        checkIfExists: true)
ch_riparian_script        = file("${projectDir}/r_scripts/riparian_paf.R",                       checkIfExists: true)
ch_assembly_report_script = file("${projectDir}/py_scripts/generate_assembly_report.py",         checkIfExists: true)
ch_coverage_book_script   = file("${projectDir}/py_scripts/bigwig_genome_book.py",               checkIfExists: true)
ch_tad_book_script        = file("${projectDir}/py_scripts/make_tad_book.py",                    checkIfExists: true)
ch_compartments_script    = file("${projectDir}/py_scripts/plot_compartments_pc1_genomewide.py", checkIfExists: true)
/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Parse sample sheet -> per-sample tuple(meta, reads)
    ch_input = parseSampleSheet(params.sample_sheet)

    // ── Ensure the NCBI taxdump once, up front. Idempotent: the module skips the download
    //    when names.dmp/nodes.dmp are already present at the target path. Needed by
    //    RESOLVE_TAXONOMY and (when enabled) the decon evidence branch.
    DOWNLOAD_TAXDUMP(
        Channel.value( file(params.diamond_taxdump_dir) ),
        Channel.value( (params.diamond?.force ?: false) as boolean )
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
                  [ name         : organismName(row),
                    kingdom      : kingdomFlag(row),
                    busco_lineage: buscoLineageFor(row) ]) }

    // Per-sample identity side-channel: (sample, {taxid, name, kingdom, busco_lineage}).
    ch_sample_identity = ch_input
        .map { meta, reads -> tuple(meta.taxid?.toString(), meta.sample, meta.species) }
        .filter { taxid, sample, override -> taxid != null }
        .combine(ch_taxonomy, by: 0)                       // many samples per taxid
        .map { taxid, sample, override, tax ->
            tuple(sample, [ taxid        : taxid,
                            name         : override ?: tax.name,
                            kingdom      : tax.kingdom,
                            busco_lineage: tax.busco_lineage ]) }

    ch_sample_identity.subscribe { sample, tax ->
        log.info "  resolved '${sample}': taxid=${tax.taxid}  name='${tax.name}'  kingdom=${tax.kingdom}  busco=${tax.busco_lineage}"
    }

    /*
    ========================================================================================
        STEP 0: Setup Decontamination Databases (if enabled) & BUSCO Database
        Runs in parallel with BAM conversion and assembly
        Only executes if decontamination is requested
    ========================================================================================
    */
    if (params.decon.run_on_contigs || params.decon.run_on_scaffolds) {
        SETUP_DECONTAM_DBS(ch_taxdump)
        
        // Store outputs for later use
        ch_gxdb_dir = SETUP_DECONTAM_DBS.out.gxdb_dir
        ch_diamond_db = SETUP_DECONTAM_DBS.out.diamond_db
    }
    
    // BUSCO databases — one download per DISTINCT lineage across all samples
    // (storeDir dedupes; the per-lineage tasks run in parallel).
    ch_busco_lineages = ch_taxonomy
        .map { taxid, tax -> tax.busco_lineage }
        .unique()
    DOWNLOAD_BUSCO_DB(ch_busco_lineages)

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
        .combine( DOWNLOAD_BUSCO_DB.out.db.map { db -> tuple(db.name, db) }, by: 0 )
        .map { lineage, taxid, db -> [ (taxid): lineage ] }
        .reduce([:]) { acc, m -> acc + m }

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
    FIND_MITO_REFERENCE(ch_mito_ref_todo)

    // Per-taxid reference for the organelle step: (taxid, ref_fasta, ref_gb)
    ch_mito_ref_by_taxid = FIND_MITO_REFERENCE.out.ref_fasta
        .join( FIND_MITO_REFERENCE.out.ref_gb )        // .join is 1:1 by taxid — one ref per species

    /*
    ========================================================================================
        STEP 1: Convert BAM to FASTQ (extract hifi_bam from the reads map)
    ========================================================================================
    */
    
    BAM_TO_FASTQ(
        ch_input.filter { meta, reads -> meta.hifi }
                .map { meta, reads -> tuple(meta, reads.hifi_bam) }
    )

    
    /*
    ========================================================================================
        STEP 2: Trim Hi-C Reads (extract the Hi-C pair from the reads map)
    ========================================================================================
    */
    TRIM_HIC(
        ch_input.filter { meta, reads -> meta.hic }
                .map { meta, reads -> tuple(meta, reads.hic_r1, reads.hic_r2) }
    )

    // Optional short-read trimming (fastp): raw shotgun -> adapter/quality-trimmed, or
    // pass-through when off. Feeds the assembly + assembly-QC path; SHORTREAD_QC stays on raw.
    ch_shortread_raw = ch_input
        .filter { meta, reads -> meta.shortread }
        .map    { meta, reads -> tuple(meta, reads.sr_r1, reads.sr_r2) }

    if (params.shortread_trim) {
        TRIM_SHORTREAD(ch_shortread_raw)
        ch_shortread_reads = TRIM_SHORTREAD.out.trimmed_reads
    } else {
        ch_shortread_reads = ch_shortread_raw
    }

    /*
    ========================================================================================
        STEP 3: Combine HiFi FASTQ with trimmed Hi-C reads
    ========================================================================================
    */
    // Full per-sample read bundle for the selector + organelle + genome-size.
    // remainder:true left-joins keep samples lacking HiFi or Hi-C (null slots), so
    // short-read-only rows flow instead of being dropped by an inner join.
    // Per-modality slots: every sample gets exactly one entry per slot — its processed
    // reads if it has that modality, else a null placeholder from ch_input (immediate).
    // Plain 1:1 joins then emit each sample as soon as ITS OWN reads are ready — no
    // waiting on other samples' BAM_TO_FASTQ / TRIM_HIC / TRIM_SHORTREAD to finish.
    ch_hifi_slot = BAM_TO_FASTQ.out.fastq
        .map { meta, fq -> [ meta.sample, fq ] }
        .mix( ch_input.filter { meta, reads -> !meta.hifi }.map { meta, reads -> [ meta.sample, null ] } )

    ch_hic_slot = TRIM_HIC.out.trimmed_reads
        .map { meta, r1, r2 -> [ meta.sample, [r1, r2] ] }
        .mix( ch_input.filter { meta, reads -> !meta.hic }.map { meta, reads -> [ meta.sample, null ] } )

    ch_sr_slot = ch_shortread_reads
        .map { meta, r1, r2 -> [ meta.sample, [r1, r2] ] }
        .mix( ch_input.filter { meta, reads -> !meta.shortread }.map { meta, reads -> [ meta.sample, null ] } )

    ch_input
        .map { meta, reads -> [ meta.sample, meta ] }
        .join( ch_hifi_slot )
        .join( ch_hic_slot )
        .join( ch_sr_slot )
        .map { sample, meta, hifi_fastq, hic_pair, sr_pair ->
            def hic_r1 = hic_pair ? hic_pair[0] : null
            def hic_r2 = hic_pair ? hic_pair[1] : null
            def sr_r1  = sr_pair  ? sr_pair[0]  : null
            def sr_r2  = sr_pair  ? sr_pair[1]  : null
            tuple(meta, hifi_fastq, hic_r1, hic_r2, sr_r1, sr_r2)
        }
        .set { ch_reads_all }
        
    // Per-sample reads for assembly QC (meryl DB + mapping): HiFi FASTQ for HiFi samples,
    // the Illumina R1+R2 pair for short-read samples. Read-source-aware QC.
    ch_reads_all
        .map { meta, hifi_fastq, hic_r1, hic_r2, sr_r1, sr_r2 ->
            meta.hifi ? tuple(meta, hifi_fastq) : tuple(meta, [sr_r1, sr_r2])
        }
        .set { ch_qc_reads }

    // Genome-size estimation (jellyfish -> GenomeScope2), concurrent with assembly.
    // Reads by assembler: HiFi for the long-read path, PE for short-read.
    ch_reads_all
        .map { meta, hifi_fastq, r1, r2, sr1, sr2 ->
            def gs_reads = (meta.assembler == 'spades') ? [ sr1, sr2 ] : [ hifi_fastq ]
            tuple(meta, gs_reads)
        }
        .set { ch_gsize_input }

    ESTIMATE_GENOME_SIZE(ch_gsize_input)


    /*
    ========================================================================================
        STEP 3b: Assemble Mitochondrial Genome (MitoHiFi)
        Runs concurrently with HIFIASM — both consume HiFi reads independently
        Depends on: BAM_TO_FASTQ + FIND_MITO_REFERENCE
    ========================================================================================
    */
    ORGANELLE_ASSEMBLY(
        ch_reads_all.map { meta, hifi_fastq, hic1, hic2, sr1, sr2 -> tuple(meta, hifi_fastq, sr1, sr2) },
        ch_mito_ref_by_taxid
    )

    /*
    ========================================================================================
        STEP 4: Assemble contigs — assembler selector (hifiasm | spades)
        hifiasm: HiFi(+Hi-C) -> hap1+hap2 (diploid) / primary (haploid)
        spades:  PE short reads -> one collapsed 'primary' assembly
    ========================================================================================
    */

    CONTIG_ASSEMBLY(ch_reads_all)

    // Fork sample-level assembly into per-hap contigs: diploid -> hap1+hap2, haploid/spades -> one primary.
    CONTIG_ASSEMBLY.out.assemblies
        .flatMap { meta, fastas ->
            def hmetas = forkHaplotypeMeta(meta)           // [hap1, hap2]  or  [primary]
            def fs = (fastas instanceof List) ? fastas.sort { it.name } : [fastas]
            [hmetas, fs].transpose().collect { hm, fa -> tuple(hm, fa) }
        }
        .set { ch_contigs }                                // per-haplotype tuple(meta, fasta)

    ch_contigs
        .branch { meta, fasta ->
            hifi:      meta.hifi
            shortread: true
        }
        .set { ch_contigs_by_type }

    /*
    ====================================================================================
        STEP 4a: Remove Mitochondrial Contigs from Nuclear Assembly
        ch_contigs is per-haplotype; attach the per-sample mitogenome (one fans out to
        both haplotypes via combine). A single FILTER_MITO_CONTIGS handles both.
    ====================================================================================
    */
    ch_contigs_by_type.hifi
        .map { meta, fasta -> [ meta.sample, meta, fasta ] }
        .combine( ORGANELLE_ASSEMBLY.out.mitogenome.map { meta, mfa -> [ meta.sample, mfa ] }, by: 0 )
        .map { sample, meta, fasta, mito_fa -> tuple(meta, fasta, mito_fa) }
        .set { ch_mito_filter_input }

    FILTER_MITO_CONTIGS(ch_mito_filter_input)

    ch_mito_filtered = FILTER_MITO_CONTIGS.out.filtered   // per-haplotype tuple(meta, fasta)

    /*
    ====================================================================================
        STEP 4b: PURGE DUPLICATES (optional)
        ch_mito_filtered is per-haplotype; attach per-sample HiFi reads (combine by sample),
        one PURGE_DUPS over both haplotypes.
    ====================================================================================
    */
    if (params.run_purge_dups) {
        ch_mito_filtered
            .map { meta, fasta -> [ meta.sample, meta, fasta ] }
            .combine( BAM_TO_FASTQ.out.fastq.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
            .map { sample, meta, fasta, hifi_reads -> tuple(meta, fasta, hifi_reads) }
            .set { ch_purge_dups_input }

        PURGE_DUPS(ch_purge_dups_input)
        ch_hifiasm_output = PURGE_DUPS.out.purged_assembly      // per-haplotype (meta, fasta)
    } else {
        ch_hifiasm_output = ch_mito_filtered                    // per-haplotype (meta, fasta)
    }

    // Short-read conditioning: REDUNDANS (reduce/scaffold/gap-close) -> optional Pilon.
    ch_contigs_by_type.shortread
        .map { meta, fasta -> [ meta.sample, meta, fasta ] }
        .combine( ch_shortread_reads.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
        .map { sample, meta, fasta, r1, r2 -> tuple(meta, fasta, r1, r2) }
        .set { ch_redundans_input }

    REDUNDANS(ch_redundans_input)

    if (params.run_pilon) {
        REDUNDANS.out.assembly
            .map { meta, fasta -> [ meta.sample, meta, fasta ] }
            .combine( ch_shortread_reads.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
            .map { sample, meta, fasta, r1, r2 -> tuple(meta, fasta, r1, r2) }
            .set { ch_pilon_input }
        PILON(ch_pilon_input)
        ch_shortread_conditioned = PILON.out.assembly
    } else {
        ch_shortread_conditioned = REDUNDANS.out.assembly
    }

    /*
    ====================================================================================
        STEP 4.5: Optional Misassembly Correction of Contig Assemblies (Inspector)
    ====================================================================================
    */
    if (params.inspector_run_on_contigs) {
        ch_hifiasm_output
            .map { meta, fasta -> [ meta.sample, meta, fasta ] }
            .combine( BAM_TO_FASTQ.out.fastq.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
            .map { sample, meta, fasta, hifi_fastq ->
                def correction_params = [
                    min_depth: params.inspector_contig_min_depth,
                    min_contig_length: params.inspector_contig_min_contig_length,
                    min_contig_length_assemblyerror: params.inspector_contig_min_contig_length_assemblyerror,
                    min_assembly_error_size: params.inspector_contig_min_assembly_error_size,
                    max_assembly_error_size: params.inspector_contig_max_assembly_error_size,
                    skip_baseerror: params.inspector_contig_skip_baseerror
                ]
                tuple(meta, fasta, hifi_fastq, "contig", correction_params)
            }
            .set { ch_correction_input }

        CORRECT_MISASSEMBLIES_CONTIG(ch_correction_input)
        ch_hifi_conditioned = CORRECT_MISASSEMBLIES_CONTIG.out.corrected   // per-hap (meta, fasta)
    } else {
        ch_hifi_conditioned = ch_hifiasm_output                           // per-hap (meta, fasta)
    }

    /*
    ====================================================================================
        STEP 5: Optional Decontamination of Contig Assemblies
        HiFi-conditioned (mito-filter / purge / correct) and short-read-conditioned
        (redundans / pilon) both flow through the SAME optional decontamination — FCS-GX is
        genome-based, so short-read assemblies are screened too (per-sample taxid via
        meta.taxid; see decontaminate_assembly.nf). Short-read gets NO Inspector correction
        (long-read-only). After decontam: HiFi continues into Hi-C scaffolding / the HiFi-only
        bypass; short-read is finished (nothing to scaffold with) → straight to FINALIZE.
    ====================================================================================
    */
    ch_assemblies_for_decontam = ch_hifi_conditioned.mix(ch_shortread_conditioned)

    if (params.decon.run_on_contigs) {
        DECONTAMINATE_ASSEMBLY_CONTIG(ch_assemblies_for_decontam, ch_gxdb_dir, "contig")
        ch_decontaminated_contigs = DECONTAMINATE_ASSEMBLY_CONTIG.out.decontaminated
    } else {
        ch_decontaminated_contigs = ch_assemblies_for_decontam
    }

    ch_individual_haplotypes = ch_decontaminated_contigs.filter { meta, fasta -> !meta.shortread }  // HiFi (+ HiFi-only)
    ch_shortread_finished    = ch_decontaminated_contigs.filter { meta, fasta ->  meta.shortread }  // → FINALIZE

    // short-read already split off, so this is HiFi-only:
    ch_hifi_only_scaffolds = ch_individual_haplotypes.filter { meta, fasta -> !meta.hic }

    // HiFi-only rows have no Hi-C, so they drop out of the Hi-C scaffolding block below.
    // Their decontaminated contigs are their "scaffolds" — rejoin at gap-filling.
    ch_hifi_only_scaffolds = ch_individual_haplotypes.filter { meta, fasta -> !meta.hic }

    /*
    ====================================================================================
        STEP 6: Map Hi-C to Assemblies (contigs or decontaminated contigs)
    ====================================================================================
    */
    // Combine each haplotype with its sample's trimmed Hi-C reads (key on meta.sample —
    // TRIM_HIC carries the sample-level meta; ch_individual_haplotypes is per-haplotype)
    ch_individual_haplotypes
        .map { meta, fasta -> [ meta.sample, meta, fasta ] }
        .combine( TRIM_HIC.out.trimmed_reads.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
        .map { sample, meta, fasta, hic_r1, hic_r2 ->
            tuple(meta, fasta, hic_r1, hic_r2, "contig")
        }
        .set { ch_hic_mapping_input }

    MAP_HIC_TO_ASSEMBLY(ch_hic_mapping_input)

    // checkpoint: contig_raw_map (BAM-level only)
    MAP_HIC_TO_ASSEMBLY.out.bam
        .map { meta, stage, bam, bai -> tuple(meta, "contig_raw_map", bam, bai) }
        .set { ch_hic_raw_bam_for_qc }

    HIC_BAM_METRICS_CONTIG(ch_hic_raw_bam_for_qc)

    // Assemblies channel for the filter / scaffold joins (meta, fasta)
    ch_individual_haplotypes
        .map { meta, fasta -> tuple(meta, fasta) }
        .set { ch_assemblies_for_qc }


    /*
    ====================================================================================
        STEP 7: Filter Hi-C BAM Files
    ====================================================================================
    */
    MAP_HIC_TO_ASSEMBLY.out.bam
        .join(ch_assemblies_for_qc)
        .map { meta, stage, bam, bai, assembly_fasta ->
            tuple(meta, stage, bam, bai, assembly_fasta)
        }
        .set { ch_bam_with_assembly }

    FILTER_HIC_BAM(ch_bam_with_assembly)

    // checkpoint: contig_filtered (pairs-level + retention); no AGP here -> []
    FILTER_HIC_BAM.out.pairs
        .join(FILTER_HIC_BAM.out.parse_stats)
        .join(FILTER_HIC_BAM.out.dedup_stats)
        .map { meta, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats ->
            tuple(meta, "contig_filtered", pairs_gz, [], parse_stats, dedup_stats)
        }
        .set { ch_hic_pairs_contig_filtered_for_qc }

    HIC_PAIRS_METRICS_CONTIG(ch_hic_pairs_contig_filtered_for_qc)
    
    /*
    ====================================================================================
        STEP 8: First Round of Scaffolding with Hi-C
    ====================================================================================
    */
    FILTER_HIC_BAM.out.bam
        .join(ch_assemblies_for_qc)
        .map { meta, stage, bam, bai, assembly_fasta ->
            tuple(meta, bam, bai, assembly_fasta, "round1", params.yahs_round1)
        }
        .set { ch_scaffolding_round1_input }

    SCAFFOLD_HIC_ROUND1(ch_scaffolding_round1_input)

    // checkpoint: scaffold_space (relabel contigs->scaffolds via round1 AGP; no remap)
    FILTER_HIC_BAM.out.pairs
        .join(SCAFFOLD_HIC_ROUND1.out.agp)
        .join(FILTER_HIC_BAM.out.parse_stats)
        .join(FILTER_HIC_BAM.out.dedup_stats)
        .map { meta, stage, pairs_gz, agp, stage2, parse_stats, stage3, dedup_stats ->
            tuple(meta, "scaffold_space", pairs_gz, agp, parse_stats, dedup_stats)
        }
        .set { ch_hic_pairs_scaffold_space_for_qc }

    HIC_PAIRS_METRICS_CONTIGSCAF(ch_hic_pairs_scaffold_space_for_qc)

    /*
    ====================================================================================
        STEP 8.5: Optional Misassembly Correction of Scaffolded Assemblies (Inspector)
    ====================================================================================
    */
    if (params.inspector_run_on_scaffolds) {
        // Combine each scaffolded haplotype with its sample's HiFi reads (key on meta.sample)
        SCAFFOLD_HIC_ROUND1.out.scaffolds
            .map { meta, scaffold -> [ meta.sample, meta, scaffold ] }
            .combine( BAM_TO_FASTQ.out.fastq.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
            .map { sample, meta, scaffold, hifi_fastq ->
                def correction_params = [
                    min_depth: params.inspector_scaffold_min_depth,
                    min_contig_length: params.inspector_scaffold_min_contig_length,
                    min_contig_length_assemblyerror: params.inspector_scaffold_min_contig_length_assemblyerror,
                    min_assembly_error_size: params.inspector_scaffold_min_assembly_error_size,
                    max_assembly_error_size: params.inspector_scaffold_max_assembly_error_size,
                    skip_baseerror: params.inspector_scaffold_skip_baseerror
                ]
                tuple(meta, scaffold, hifi_fastq, "scaffold", correction_params)
            }
            .set { ch_scaffold_correction_input }

        CORRECT_MISASSEMBLIES_SCAFFOLD(ch_scaffold_correction_input)
        ch_scaffolds_for_decontam = CORRECT_MISASSEMBLIES_SCAFFOLD.out.corrected   // per-hap (meta, fasta)
    } else {
        ch_scaffolds_for_decontam = SCAFFOLD_HIC_ROUND1.out.scaffolds              // per-hap (meta, fasta)
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
        ================================================================================
            STEP 10: Map Hi-C to Final Scaffolds
        ================================================================================
        */
        ch_final_scaffolds
            .map { meta, fasta -> [ meta.sample, meta, fasta ] }
            .combine( TRIM_HIC.out.trimmed_reads.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
            .map { sample, meta, fasta, hic_r1, hic_r2 ->
                tuple(meta, fasta, hic_r1, hic_r2, "scaffold")
            }
            .set { ch_hic_scaffold_mapping_input }

        MAP_HIC_TO_SCAFFOLD(ch_hic_scaffold_mapping_input)

        // checkpoint: scaffold_round2_raw_map (BAM-level only)
        MAP_HIC_TO_SCAFFOLD.out.bam
            .map { meta, stage, bam, bai -> tuple(meta, "scaffold_round2_raw_map", bam, bai) }
            .set { ch_hic_scaffold_raw_bam_for_qc }

        HIC_BAM_METRICS_SCAFFOLD(ch_hic_scaffold_raw_bam_for_qc)

        /*
        ================================================================================
            STEP 11: Filter Hi-C BAM mapped to final scaffolds
        ================================================================================
        */
        MAP_HIC_TO_SCAFFOLD.out.bam
            .join(ch_final_scaffolds)
            .map { meta, stage, bam, bai, scaffold_fasta ->
                tuple(meta, stage, bam, bai, scaffold_fasta)
            }
            .set { ch_bam_with_scaffold }

        FILTER_HIC_BAM_SCAFFOLD(ch_bam_with_scaffold)

        // checkpoint: scaffold_round2_filtered (pairs-level + retention); already in scaffold1 names -> []
        FILTER_HIC_BAM_SCAFFOLD.out.pairs
            .join(FILTER_HIC_BAM_SCAFFOLD.out.parse_stats)
            .join(FILTER_HIC_BAM_SCAFFOLD.out.dedup_stats)
            .map { meta, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats ->
                tuple(meta, "scaffold_round2_filtered", pairs_gz, [], parse_stats, dedup_stats)
            }
            .set { ch_hic_pairs_scaffold_round2_filtered_for_qc }

        HIC_PAIRS_METRICS_SCAFFOLD(ch_hic_pairs_scaffold_round2_filtered_for_qc)


        /*
        ================================================================================
            STEP 12: Second Round of Scaffolding
        ================================================================================
        */
        FILTER_HIC_BAM_SCAFFOLD.out.bam
            .join(ch_final_scaffolds)
            .map { meta, stage, bam, bai, scaffold_fasta ->
                tuple(meta, bam, bai, scaffold_fasta, "round2", params.yahs_round2)
            }
            .set { ch_second_scaffolding_input }

        SCAFFOLD_HIC_ROUND2(ch_second_scaffolding_input)
        ch_final_scaffolds_round2 = SCAFFOLD_HIC_ROUND2.out.scaffolds   // per-hap (meta, fasta)

        // checkpoint: scaffold_round2_space (relabel scaffold1->scaffold2 via round2 AGP; no remap)
        FILTER_HIC_BAM_SCAFFOLD.out.pairs
            .join(SCAFFOLD_HIC_ROUND2.out.agp)
            .join(FILTER_HIC_BAM_SCAFFOLD.out.parse_stats)
            .join(FILTER_HIC_BAM_SCAFFOLD.out.dedup_stats)
            .map { meta, stage, pairs_gz, agp, stage2, parse_stats, stage3, dedup_stats ->
                tuple(meta, "scaffold_round2_space", pairs_gz, agp, parse_stats, dedup_stats)
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

    // HiFi-only assemblies are NOT gap-filled — no Hi-C scaffolding means no scaffold gaps to
    // close. They rejoin the finishing chain at teloclip/finalize below (they still have HiFi
    // reads), mirroring how short-read rejoins at ch_final_assembly.
    
    // Combine scaffolds with sample HiFi reads for gap filling (key on meta.sample)
    ch_scaffolds_for_gap_filling
        .map { meta, scaffold -> [ meta.sample, meta, scaffold ] }
        .combine( BAM_TO_FASTQ.out.fastq.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
        .map { sample, meta, scaffold, hifi_fastq -> tuple(meta, scaffold, hifi_fastq) }
        .set { ch_gap_filling_input }
    
    // Run gap filling
    GAP_FILLING(ch_gap_filling_input)

    // Gap-filled Hi-C scaffolds + HiFi-only contigs (which correctly skipped gap-fill) both
    // continue to teloclip/finalize.
    ch_post_gap_fill = GAP_FILLING.out.filled_assembly.mix(ch_hifi_only_scaffolds)

    /*
    ========================================================================================
        STEP 13b: Teloclip — Extend scaffolds with missing telomeres (Optional)
        Maps raw HiFi reads back to gap-filled scaffolds to find soft-clipped
        alignments at scaffold ends containing telomeric motifs, then appends
        the overhang sequence to recover missing telomeres.
    ========================================================================================
    */
    if (params.run_teloclip_extend) {
        // Combine gap-filled assemblies with sample HiFi reads (key on meta.sample)
        ch_post_gap_fill
            .map { meta, filled_fa -> [ meta.sample, meta, filled_fa ] }
            .combine( BAM_TO_FASTQ.out.fastq.map { meta, fq -> [ meta.sample, fq ] }, by: 0 )
            .map { sample, meta, filled_fa, hifi_fastq -> tuple(meta, filled_fa, hifi_fastq) }
            .set { ch_teloclip_input }

        TELOCLIP_EXTEND(ch_teloclip_input)

        // Collect teloclip stats across all haplotypes
        COLLECT_TELOCLIP_STATS(
            TELOCLIP_EXTEND.out.stats.map { meta, stats -> stats }.collect()
        )

        // The teloclip-extended assembly becomes the "final" assembly
        ch_final_assembly = TELOCLIP_EXTEND.out.extended_assembly
        ch_teloclip_stats_for_report = COLLECT_TELOCLIP_STATS.out.stats.ifEmpty(file('NO_TELOCLIP'))
    } else {
        // No teloclip — gap-filled assembly IS the final assembly
        ch_final_assembly = ch_post_gap_fill
        ch_teloclip_stats_for_report = Channel.of(file('NO_TELOCLIP'))
    }

    // =========================================================================
    //  FINALIZE ASSEMBLY — now uses ch_final_assembly (post-teloclip if enabled)
    // =========================================================================
    ch_final_assembly = ch_final_assembly.mix(ch_shortread_finished)

    FINALIZE_ASSEMBLY(ch_final_assembly)
    ch_finalized_assembly = FINALIZE_ASSEMBLY.out.assembly

    /*
    ========================================================================================
        STEP 14: Finalization
    ========================================================================================
    */
    // 1. Contact Maps for Final Assemblies
    //    REPLACE: GAP_FILLING.out.filled_assembly → ch_final_assembly
    if (params.make_final_contact_maps) {
        // Combine final per-hap assemblies with sample Hi-C reads (key on meta.sample)
        ch_finalized_assembly
            .map { meta, final_fa -> [ meta.sample, meta, final_fa ] }
            .combine( TRIM_HIC.out.trimmed_reads.map { meta, r1, r2 -> [ meta.sample, r1, r2 ] }, by: 0 )
            .map { sample, meta, final_fa, hic_r1, hic_r2 ->
                tuple(meta, final_fa, hic_r1, hic_r2, 'final')
            }
            .set { ch_final_hic_map_input }

        MAP_HIC_TO_FINAL(ch_final_hic_map_input)

        // checkpoint: final_raw_map (BAM-level only)
        MAP_HIC_TO_FINAL.out.bam
            .map { meta, stage, bam, bai -> tuple(meta, "final_raw_map", bam, bai) }
            .set { ch_final_raw_bam_qc }

        HIC_BAM_METRICS_FINAL(ch_final_raw_bam_qc)

        MAP_HIC_TO_FINAL.out.bam
            .join(ch_finalized_assembly)
            .map { meta, stage, bam, bai, final_fa ->
                tuple(meta, "final", bam, bai, final_fa)
            }
            .set { ch_final_filter_input }

        FILTER_HIC_BAM_FINAL(ch_final_filter_input)

        // checkpoint: final_filtered (pairs-level + retention); already in final names -> []
        FILTER_HIC_BAM_FINAL.out.pairs
            .join(FILTER_HIC_BAM_FINAL.out.parse_stats)
            .join(FILTER_HIC_BAM_FINAL.out.dedup_stats)
            .map { meta, stage, pairs_gz, stage2, parse_stats, stage3, dedup_stats ->
                tuple(meta, "final_filtered", pairs_gz, [], parse_stats, dedup_stats)
            }
            .set { ch_final_pairs_qc }

        HIC_PAIRS_METRICS_FINAL(ch_final_pairs_qc)

        FILTER_HIC_BAM_FINAL.out.pairs
            .join(ch_finalized_assembly)
            .map { meta, stage, pairs_gz, final_fasta ->
                tuple(meta, pairs_gz, final_fasta, "final")
            }
            .set { ch_contact_map_final_input }

        CONTACT_MAP_FINAL(ch_contact_map_final_input)

        if (params.hic_balance) {
            HIC_COMPARTMENTS(
                CONTACT_MAP_FINAL.out.mcool,
                params.compartment_resolution ?: 250000,
                params.compartment_min_contig_bp ?: 5000000,
                params.compartment_max_contigs ?: 30,
                ch_compartments_script
            )

            HIC_TADS(
                CONTACT_MAP_FINAL.out.mcool,
                params.tad_resolution ?: 50000,
                params.tad_window_bp ?: 500000,
                params.tad_min_contig_bp ?: 5000000,
                params.tad_max_contigs ?: 0,
                ch_tad_book_script
            )
        }
    }

    // String-id view of the finalized per-hap assemblies for the leaf viz/QUAST steps.
    // These only need an id (for pairing, labels, output naming) — meta itself isn't used.
    ch_finalized_assembly
        .map { meta, fa -> tuple(meta.id, fa) }
        .set { ch_final_by_id }

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
        ch_final_by_id 
            .toSortedList { a, b -> a[0] <=> b[0] }
            .flatMap { assemblies ->
                def pairs = []

                if (params.pairwise_alignment_mode == 'within_sample') {
                    // Only compare hap1 vs hap2 within each sample.
                    // A haploid sample has a single 'primary' assembly and therefore no
                    // within-sample comparison — it is intentionally excluded below by the
                    // haps.size() == 2 check (its group has size 1).
                    def grouped = assemblies.groupBy { haplotype_id, fasta ->
                        haplotype_id.replaceAll(/_(hap[12]|primary)$/, '')
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
        PAIRWISE_ALIGNMENT(ch_pairwise_input, SETUP_PAFR.out.ready, ch_dotplot_script)

        // Riparian plot input — uses ch_final_assembly
        ch_paf_with_asm1 = PAIRWISE_ALIGNMENT.out.paf
            .map { hap1, hap2, paf -> tuple(hap1, hap2, paf) }
            .combine(ch_final_by_id, by: 0)

        ch_riparian_input = ch_paf_with_asm1
            .map { hap1, hap2, paf, fasta1 -> tuple(hap2, hap1, fasta1, paf) }
            .combine(ch_final_by_id, by: 0)
            .map { hap2, hap1, fasta1, paf, fasta2 ->
                tuple(hap1, fasta1, hap2, fasta2, paf)
            }

        RIPARIAN_PLOT(ch_riparian_input, ch_riparian_script)

        COLLECT_PAIRWISE_RESULTS(
            PAIRWISE_ALIGNMENT.out.qc.map { id1, id2, qc_file -> qc_file }.collect()
        )
        ch_pairwise_summary = COLLECT_PAIRWISE_RESULTS.out.summary.ifEmpty(file('NO_PAIRWISE'))
    } else {
        ch_pairwise_summary = Channel.of(file('NO_PAIRWISE'))
    }

    /*
    ========================================================================================
        QUAST Final Comparison - All Gap-Filled Genomes
    ========================================================================================
    */
    // Collect all gap-filled assemblies with their labels for cross-sample comparison
    ch_final_by_id
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
    TIDK_EXPLORE(ch_final_by_id)

    // 2. Search: windowed telomere repeat quantification
    TIDK_SEARCH(ch_final_by_id)

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
    BUILD_MERYL_DB(ch_qc_reads)

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
        ch_input.filter { meta, reads -> meta.hic }
                .map { meta, reads -> tuple(meta, reads.hic_r1, reads.hic_r2) },
        "raw"
    )

    /*
    ========================================================================================
        QC HiFi Reads
    ========================================================================================
    */
    HIFI_QC(
        BAM_TO_FASTQ.out.fastq
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

    // Short-read input QC — raw + trimmed (mirrors HIC_QC_RAW / HIC_QC_TRIMMED)
    SHORTREAD_QC_RAW(
        ch_input.filter { meta, reads -> meta.shortread }
                .map { meta, reads -> tuple(meta, reads.sr_r1, reads.sr_r2) },
        "raw"
    )
    if (params.shortread_trim) {
        SHORTREAD_QC_TRIMMED(TRIM_SHORTREAD.out.trimmed_reads, "trimmed")
    }
    
    /*
    ========================================================================================
        Assembly QC
    ========================================================================================
    */

    // QC stage selector: 'all_stages' (default) runs QC at every checkpoint;
    // 'final_only' runs QC on the final assembly only.
    def run_all_qc = ( (params.qc_mode ?: 'all_stages') != 'final_only' )

     // QC raw hifiasm contigs (per-hap fork — was HIFIASM.out.assemblies triple)
    if (run_all_qc) {
        ASSEMBLY_QC_INITIAL(
            ch_contigs,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'contig'
        )
    }

    // QC mito-filtered contigs
    if (run_all_qc) {
        ASSEMBLY_QC_MITO_FILTERED(
            ch_mito_filtered,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'contig_mito_filtered'
        )
    }

    // QC purged contigs
    if (run_all_qc && params.run_purge_dups) {
        ASSEMBLY_QC_PURGED(
            ch_hifiasm_output,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'contig_purged'
        )
    }

    // QC redundans-reduced short-read contigs — short-read analogue of purge_dups, reported
    // at the same 'contig_purged' (ctg.purged) stage. No-op when there are no short-read samples.
    if (run_all_qc) {
        ASSEMBLY_QC_REDUNDANS(
            ch_shortread_conditioned, ch_qc_reads, BUILD_MERYL_DB.out.meryl_db, ch_busco_db,
            'contig_purged'
        )
    }

    // QC Inspector-corrected contigs
    if (run_all_qc && params.inspector_run_on_contigs) {
        ASSEMBLY_QC_CONTIG_CORRECTED(
            CORRECT_MISASSEMBLIES_CONTIG.out.corrected,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'contig_corrected'
        )
    }

    // QC decontaminated contigs
    if (run_all_qc && params.decon.run_on_contigs) {
        ASSEMBLY_QC_CONTIG_DECONTAM(
            DECONTAMINATE_ASSEMBLY_CONTIG.out.decontaminated,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'contig_decontam'
        )
    }

    // QC round-1 scaffolds
    if (run_all_qc) {
        ASSEMBLY_QC_SCAFFOLD(
            SCAFFOLD_HIC_ROUND1.out.scaffolds,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'scaffold'
        )
    }

    // QC Inspector-corrected scaffolds
    if (run_all_qc && params.inspector_run_on_scaffolds) {
        ASSEMBLY_QC_SCAFFOLD_CORRECTED(
            CORRECT_MISASSEMBLIES_SCAFFOLD.out.corrected,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'scaffold_corrected'
        )
    }

    // QC decontaminated scaffolds
    if (run_all_qc && params.decon.run_on_scaffolds) {
        ASSEMBLY_QC_SCAFFOLD_DECONTAM(
            DECONTAMINATE_ASSEMBLY_SCAFFOLD.out.decontaminated,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'scaffold_decontam'
        )
    }

    // QC round-2 scaffolds
    if (run_all_qc && params.run_scaffold_round2) {
        ASSEMBLY_QC_SCAFFOLD_ROUND2(
            ch_final_scaffolds_round2,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'scaffold_round2'
        )
    }

    // QC gap-filled genomes — intermediate HiFi-path stage (short-read has no gap-fill).
    if (run_all_qc) {
        ASSEMBLY_QC_GAP_FILLED(
            GAP_FILLING.out.filled_assembly,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'gap_filled'
        )
    }

    // QC teloclip-extended assembly (pre-FINALIZE) — kept independent of the final QC so the
    // effect of FINALIZE (scaffold renaming, future small-contig trimming, etc.) stays
    // visible. HiFi-path only, when teloclip is enabled.
    if (run_all_qc && params.run_teloclip_extend) {
        ASSEMBLY_QC_TELOCLIP(
            TELOCLIP_EXTEND.out.extended_assembly,
            ch_qc_reads,
            BUILD_MERYL_DB.out.meryl_db,
            ch_busco_db,
            'teloclip_extended'
        )
    }

    // Final QC — ALWAYS runs on the finalized assembly (post-FINALIZE), independent of
    // teloclip and of qc_mode (final_only needs it too). The 'final' stage for EVERY sample.
    ASSEMBLY_QC_FINAL(
        ch_finalized_assembly,
        ch_qc_reads,
        BUILD_MERYL_DB.out.meryl_db,
        ch_busco_db,
        'final'
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
                ch_qc_reads,
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
                ch_qc_reads,
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
    
    // Collect all assembly QC summaries (respecting qc_mode + which steps ran)
    ch_all_assembly_summaries = Channel.empty()

    if (run_all_qc) {
        ch_all_assembly_summaries = ch_all_assembly_summaries
            .mix(ASSEMBLY_QC_INITIAL.out.assembly_summary)
            .mix(ASSEMBLY_QC_MITO_FILTERED.out.assembly_summary)
            .mix(ASSEMBLY_QC_SCAFFOLD.out.assembly_summary)
            .mix(ASSEMBLY_QC_REDUNDANS.out.assembly_summary)

        if (params.run_purge_dups)
            ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_PURGED.out.assembly_summary)
        if (params.inspector_run_on_contigs)
            ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_CONTIG_CORRECTED.out.assembly_summary)
        if (params.decon.run_on_contigs)
            ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_CONTIG_DECONTAM.out.assembly_summary)
        if (params.inspector_run_on_scaffolds)
            ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_SCAFFOLD_CORRECTED.out.assembly_summary)
        if (params.decon.run_on_scaffolds)
            ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_SCAFFOLD_DECONTAM.out.assembly_summary)
        if (params.run_scaffold_round2)
            ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_SCAFFOLD_ROUND2.out.assembly_summary)
    }

    // Gap-filled summary — intermediate HiFi-path stage.
    if (run_all_qc) {
        ch_all_assembly_summaries = ch_all_assembly_summaries
            .mix(ASSEMBLY_QC_GAP_FILLED.out.assembly_summary)
    }

    // Teloclip-extended summary — intermediate; only when teloclip ran under run_all_qc.
    if (run_all_qc && params.run_teloclip_extend) {
        ch_all_assembly_summaries = ch_all_assembly_summaries
            .mix(ASSEMBLY_QC_TELOCLIP.out.assembly_summary)
    }

    // Final summary — always (ASSEMBLY_QC_FINAL runs unconditionally).
    ch_all_assembly_summaries = ch_all_assembly_summaries
        .mix(ASSEMBLY_QC_FINAL.out.assembly_summary)

    ch_final_busco = ASSEMBLY_QC_FINAL.out.busco_results

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
        ch_all_assembly_summaries.map { sample_id, qc_label, tsv -> tsv }.collect().ifEmpty([]),
        ch_all_bam_metrics.map { meta, checkpoint, tsv -> tsv }.collect().ifEmpty([]),
        ch_all_pairs_metrics.map { meta, checkpoint, tsv -> tsv }.collect().ifEmpty([]),
        ch_compile_qc_script
    )

    // Make interactive HTML assembly viewer
    ASSEMBLY_REPORT(COMPILE_FINAL_QC.out.metrics, ch_assembly_report_script)

    /*
    ========================================================================================
    SNAIL PLOTS FOR FINAL ASSEMBLIES
    ========================================================================================
    */
    // Join gap-filled assemblies with their BUSCO results
    // BUSCO output from ASSEMBLY_QC_TELOCLIP is per-haplotype
    ch_finalized_assembly
        .join(ch_final_busco)
        .map { meta, assembly, busco_dir ->
            tuple(meta.id, assembly, busco_dir, "final")
        }
        .set { ch_snail_plot_final_input }

    SNAIL_PLOT_FINAL(ch_snail_plot_final_input)
    /*
    ========================================================================================
        COVERAGE BOOK - HiFi coverage visualization for final assemblies
    ========================================================================================
    */

    // Combine gap-filled assemblies with HiFi reads for coverage book
    ch_finalized_assembly
        .map { meta, final_fa -> tuple(meta.sample, meta, final_fa) }
        .combine(
            BAM_TO_FASTQ.out.fastq.map { meta, hifi_fastq -> tuple(meta.sample, hifi_fastq) },
            by: 0
        )
        .map { sample, meta, final_fa, hifi_fastq ->
            tuple(meta, final_fa, hifi_fastq)
        }
        .set { ch_coverage_book_input }

    COVERAGE_BOOK(ch_coverage_book_input, ch_coverage_book_script)

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
        .map { meta, fasta -> "assembly\t${meta.id}\t.\t${fasta.name}\tassembly/final" }

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
            .flatMap { meta, stage, pngs ->
                def png_list = pngs instanceof List ? pngs : [pngs]
                png_list.collect { png ->
                    "contact_map\t${meta.id}\t.\t${png.name}\tcontact_maps"
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
    ch_manifest_mito_gb = ORGANELLE_ASSEMBLY.out.annotation
        .map { meta, gb -> "mito_genbank\t${meta.sample}\t.\t${gb.name}\tmitogenome" }

    ch_manifest_mito_stats = ORGANELLE_ASSEMBLY.out.stats
        .map { meta, tsv -> "mito_stats\t${meta.sample}\t.\t${tsv.name}\tmitogenome" }

    ch_manifest_mito_circular = ORGANELLE_ASSEMBLY.out.circular_map
        .map { meta, png -> "mito_gene_map\t${meta.sample}\t.\t${png.name}\tmitogenome" }

    // ---- Combine all manifest entries into a single TSV ----
    ch_manifest_assemblies
        .mix(ch_manifest_snails)
        .mix(ch_manifest_contact_maps)
        .mix(ch_manifest_dotplots)
        .mix(ch_manifest_riparian)
        .mix(ch_manifest_compiled_csv)
        .mix(ch_manifest_qc_plots)
        .mix(ch_manifest_assembly_report)
        .mix(ch_manifest_mito_gb)
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
    ch_mito_stats_for_report = ORGANELLE_ASSEMBLY.out.stats
        .map { meta, tsv -> tsv }
        .collectFile(
            name: 'all_mito_stats.tsv',
            keepHeader: true,
            skip: 1
        )
        .ifEmpty(file('NO_MITO_STATS'))

    // ---- Per-sample taxonomy + genome-size estimate (4b-i Increment 4) ----
    ch_sample_taxonomy_tsv = ch_sample_identity
        .map { sample, tax ->
            "${sample}\t${tax.taxid}\t${tax.name}\t${tax.kingdom}\t${tax.busco_lineage}"
        }
        .collectFile(name: 'sample_taxonomy.tsv',
                     seed: 'sample\ttaxid\tspecies\tkingdom\tbusco_lineage',
                     newLine: true)
        .ifEmpty(file('NO_TAXONOMY'))

    ch_genome_size_tsv = ESTIMATE_GENOME_SIZE.out.size
        .map { meta, size_file -> "${meta.sample}\t${size_file.text.trim()}" }
        .collectFile(name: 'genome_sizes.tsv',
                     seed: 'sample\test_genome_size_bp',
                     newLine: true)
        .ifEmpty(file('NO_GENOME_SIZE'))

    // ---- Call SUMMARY_REPORT ----
    SUMMARY_REPORT(
        ch_report_manifest,
        COMPILE_FINAL_QC.out.metrics,
        ch_telomere_for_report,
        ch_pairwise_summary,
        ch_mito_stats_for_report,
        ch_teloclip_stats_for_report,
        ch_sample_taxonomy_tsv,
        ch_genome_size_tsv,
        ch_summary_report_script       
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