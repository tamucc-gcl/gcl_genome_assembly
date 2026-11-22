nextflow.enable.dsl = 2

include { CONTAMINATION_REBUILD_REPORT; CONTAMINATION_DIAGNOSTICS; CONTAMINATION_BLOBDIR
} from '../../modules/contamination_evidence.nf'
include { CONTAMINATION_SEARCH } from '../../workflows/contamination_search.nf'
include { DOWNLOAD_CONTAMINATION_NT; VERIFY_CONTAMINATION_NT } from '../../modules/download_contamination_nt.nf'

params.fasta = null
params.assembly = null
params.stage = 'post_fcs'
params.taxid = null
params.host_rank = 'phylum'
params.taxonomy_dir = null
params.protein_db_dir = null
params.nucleotide_db = null
params.nt_db_root = null
params.nt_release = null
params.nt_collection = 'nt'
params.targeted_blastn = false
params.contamination_followup_max_bp = 100000
params.coverage_bam = null
params.coverage_tsv = null
params.previous_report = null
params.window_bp = 100000
params.chunk_bp = 5000000
params.search_mode = 'full' // Legacy standalone default; use staged for bounded post-FCS pilots.
params.contamination_sample_window_bp = 10000
params.contamination_sample_spacing_bp = 1000000
params.contamination_sample_short_contig_bp = 100000
params.contamination_sample_gc_delta = 0.15
params.contamination_sample_depth_fold = 4
params.contamination_blastx_max_forks = 8
params.contamination_blastn_max_forks = 8
params.contamination_megablast_execution = 'chunks'
params.contamination_megablast_nodes = 7
params.contamination_megablast_slots = 8
params.contamination_megablast_threads = 8
params.contamination_megablast_scratch = '/scratch'
params.contamination_megablast_time = '96h'
params.outdir = 'contamination-evidence'

workflow {
    if (!params.fasta || !params.taxonomy_dir)
        error 'Supply --fasta and --taxonomy_dir'
    def fasta = file(params.fasta, checkIfExists: true)
    def taxdump = file(params.taxonomy_dir, checkIfExists: true).toAbsolutePath().toString()
    def helpers = ['contamination_review.py', 'contamination_evidence.py'].collect {
        file("${projectDir}/../../py_scripts/${it}", checkIfExists: true)
    }
    def plotter = file("${projectDir}/../../r_scripts/plot_contamination_diagnostics.R", checkIfExists: true)

    if (params.previous_report) {
        def previous = file(params.previous_report, checkIfExists: true)
        def metadata = new groovy.json.JsonSlurper().parseText(previous.resolve('evidence.json').text)
        CONTAMINATION_REBUILD_REPORT(previous, fasta, helpers, metadata.assembly)
        CONTAMINATION_DIAGNOSTICS(CONTAMINATION_REBUILD_REPORT.out, plotter, false)
        CONTAMINATION_BLOBDIR(CONTAMINATION_REBUILD_REPORT.out.map { id, review ->
            tuple(id, fasta, review, metadata.coverage != 'unavailable', metadata.host_taxid)
        }, taxdump, helpers, false)
    } else {
        if (!(params.search_mode in ['staged', 'full']))
            error 'Set --search_mode to staged or full'
        if (params.search_mode == 'staged' && params.stage != 'post_fcs')
            error 'Staged review requires --stage post_fcs'
        ['contamination_blastx_max_forks', 'contamination_blastn_max_forks', 'contamination_followup_max_bp'].each { name ->
            def value = params[name].toString()
            if (!value.isInteger() || value.toInteger() < 1)
                error "Supply a positive integer for --${name}"
        }
        if (!params.assembly || !params.taxid || !params.protein_db_dir)
            error 'Supply --assembly, --taxid and --protein_db_dir for a new search'
        if (!(params.taxid.toString() ==~ /[1-9][0-9]*/) ||
            !(params.window_bp.toString() ==~ /[1-9][0-9]*/) ||
            !(params.chunk_bp.toString() ==~ /[1-9][0-9]*/) ||
            params.chunk_bp.toLong() < params.window_bp.toLong())
            error 'Supply a positive taxid and window size, with chunk size at least as large as the window'
        if (params.coverage_bam && params.coverage_tsv)
            error 'Supply either --coverage_bam or --coverage_tsv'
        def managed = !!(params.nt_db_root || params.nt_release)
        if (managed && (params.nucleotide_db || !params.nt_db_root ||
            !(params.nt_release?.toString() ==~ /\d{4}-\d{2}-\d{2}/) ||
            !(params.nt_collection.toString() ==~ /[A-Za-z0-9_]+/)))
            error 'Supply --nucleotide_db, or --nt_db_root and --nt_release YYYY-MM-DD for a managed NCBI database'
        def ntDb = params.nucleotide_db ? file(params.nucleotide_db).toAbsolutePath().toString() : ''
        def ntManifest = []
        if (managed) {
            def root = file(params.nt_db_root).toAbsolutePath().toString()
            def downloader = file("${projectDir}/../../py_scripts/download_contamination_nt.py", checkIfExists: true)
            DOWNLOAD_CONTAMINATION_NT(root, params.nt_collection, params.nt_release, downloader)
            VERIFY_CONTAMINATION_NT(DOWNLOAD_CONTAMINATION_NT.out.bundle, root,
                params.nt_collection, params.nt_release, downloader)
            ntDb = VERIFY_CONTAMINATION_NT.out.prefix.collect().map { prefixes -> prefixes[0] }
            ntManifest = VERIFY_CONTAMINATION_NT.out.manifest
        }
        def mega = managed || !!params.nucleotide_db
        def targeted = params.targeted_blastn.toString().toBoolean()
        if (targeted && !mega)
            error 'Targeted BLASTN requires a nucleotide database'
        def bam = params.coverage_bam ? file(params.coverage_bam, checkIfExists: true) : []
        def coverage = params.coverage_tsv ? file(params.coverage_tsv, checkIfExists: true) : []
        CONTAMINATION_SEARCH(
            Channel.of(tuple(params.assembly, fasta, params.taxid)),
            Channel.of(tuple(params.assembly, bam, [], coverage)),
            taxdump, [protein: params.protein_db_dir, stage: params.stage, window: params.window_bp,
                chunk: params.chunk_bp, mode: params.search_mode, nt: ntDb, mega: mega, targeted: targeted,
                followup: params.contamination_followup_max_bp, rank: params.host_rank, nested: false],
            file("${projectDir}/../.."), ntManifest)
    }
}
