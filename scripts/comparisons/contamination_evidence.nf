nextflow.enable.dsl = 2

include { CONTAMINATION_PREPARE; CONTAMINATION_BLASTX; CONTAMINATION_SELECT_NUCLEOTIDE;
          CONTAMINATION_BLASTN as CONTAMINATION_MEGABLAST;
          CONTAMINATION_BLASTN as CONTAMINATION_TARGETED_BLASTN;
          CONTAMINATION_REPORT; CONTAMINATION_DIAGNOSTICS;
          CONTAMINATION_BLOBDIR } from '../../modules/contamination_evidence.nf'
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
params.targeted_blastn = true
params.coverage_bam = null
params.coverage_tsv = null
params.window_bp = 100000
params.chunk_bp = 5000000
params.outdir = 'contamination-evidence'

workflow {
    if (!params.fasta || !params.assembly || !params.taxid ||
        !params.taxonomy_dir || !params.protein_db_dir)
        error 'Supply --fasta, --assembly, --taxid, --taxonomy_dir and --protein_db_dir'
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
    def fasta = file(params.fasta, checkIfExists: true)
    def bundle = file(params.protein_db_dir, checkIfExists: true).toAbsolutePath()
    def database = file(bundle.resolve('reference_proteomes.dmnd'), checkIfExists: true).toString()
    def manifest = file(bundle.resolve('protein_database.json'), checkIfExists: true)
    def taxdump = file(params.taxonomy_dir, checkIfExists: true).toAbsolutePath().toString()
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
    def targeted = mega && params.targeted_blastn.toString().toBoolean()
    def bam = params.coverage_bam ? file(params.coverage_bam, checkIfExists: true) : []
    def coverage = params.coverage_tsv ? file(params.coverage_tsv, checkIfExists: true) : []
    def helpers = ['contamination_review.py', 'contamination_evidence.py'].collect {
        file("${projectDir}/../../py_scripts/${it}", checkIfExists: true)
    }
    def plotter = file("${projectDir}/../../r_scripts/plot_contamination_diagnostics.R", checkIfExists: true)

    CONTAMINATION_PREPARE(fasta, params.assembly, params.stage, params.window_bp,
        params.chunk_bp, database, manifest, taxdump, helpers)
    def chunks = CONTAMINATION_PREPARE.out.chunks.flatten()
    CONTAMINATION_BLASTX(chunks, database)
    def hits = CONTAMINATION_BLASTX.out.collect().ifEmpty([])
    def megaHits = []
    if (mega) {
        CONTAMINATION_MEGABLAST(chunks, ntDb, 'megablast')
        megaHits = CONTAMINATION_MEGABLAST.out.collect().ifEmpty([])
    }
    CONTAMINATION_SELECT_NUCLEOTIDE(fasta, CONTAMINATION_PREPARE.out.queries, hits, megaHits, mega,
        taxdump, params.taxid, params.host_rank, helpers)
    def nucleotide = []
    if (targeted) {
        CONTAMINATION_TARGETED_BLASTN(CONTAMINATION_SELECT_NUCLEOTIDE.out.chunks.flatten(), ntDb, 'blastn')
        nucleotide = CONTAMINATION_TARGETED_BLASTN.out.collect().ifEmpty([])
    }
    CONTAMINATION_REPORT(CONTAMINATION_PREPARE.out.queries, CONTAMINATION_PREPARE.out.metrics,
        CONTAMINATION_SELECT_NUCLEOTIDE.out.selected, nucleotide, ntDb, targeted, ntManifest,
        bam, coverage, taxdump,
        params.taxid, params.host_rank, manifest, helpers)
    CONTAMINATION_DIAGNOSTICS(CONTAMINATION_REPORT.out, plotter)
    CONTAMINATION_BLOBDIR(fasta, CONTAMINATION_REPORT.out,
        !!(params.coverage_bam || params.coverage_tsv), taxdump, params.taxid, helpers)
}
