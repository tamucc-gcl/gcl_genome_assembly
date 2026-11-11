nextflow.enable.dsl = 2

include { CONTAMINATION_PREPARE; CONTAMINATION_BLASTX; CONTAMINATION_SELECT_NUCLEOTIDE;
          CONTAMINATION_BLASTN; CONTAMINATION_REPORT; CONTAMINATION_DIAGNOSTICS;
          CONTAMINATION_BLOBDIR } from '../../modules/contamination_evidence.nf'

params.fasta = null
params.assembly = null
params.stage = 'post_fcs'
params.taxid = null
params.host_rank = 'phylum'
params.taxonomy_dir = null
params.protein_db_dir = null
params.nucleotide_db = null
params.coverage_bam = null
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
    def fasta = file(params.fasta, checkIfExists: true)
    def bundle = file(params.protein_db_dir, checkIfExists: true).toAbsolutePath()
    def database = file(bundle.resolve('reference_proteomes.dmnd'), checkIfExists: true).toString()
    def manifest = file(bundle.resolve('protein_database.json'), checkIfExists: true)
    def taxdump = file(params.taxonomy_dir, checkIfExists: true).toAbsolutePath().toString()
    def ntDb = params.nucleotide_db ? file(params.nucleotide_db).toAbsolutePath().toString() : ''
    def bam = params.coverage_bam ? file(params.coverage_bam, checkIfExists: true) : []
    def helpers = ['contamination_review.py', 'contamination_evidence.py'].collect {
        file("${projectDir}/../../py_scripts/${it}", checkIfExists: true)
    }
    def plotter = file("${projectDir}/../../r_scripts/plot_contamination_diagnostics.R", checkIfExists: true)

    CONTAMINATION_PREPARE(fasta, params.assembly, params.stage, params.window_bp,
        params.chunk_bp, database, manifest, taxdump, helpers)
    CONTAMINATION_BLASTX(CONTAMINATION_PREPARE.out.chunks.flatten(), database)
    def hits = CONTAMINATION_BLASTX.out.collect().ifEmpty([])
    CONTAMINATION_SELECT_NUCLEOTIDE(fasta, CONTAMINATION_PREPARE.out.queries, hits,
        taxdump, params.taxid, params.host_rank, helpers)
    def nucleotide = []
    if (ntDb) {
        CONTAMINATION_BLASTN(CONTAMINATION_SELECT_NUCLEOTIDE.out, ntDb)
        nucleotide = CONTAMINATION_BLASTN.out
    }
    CONTAMINATION_REPORT(CONTAMINATION_PREPARE.out.queries, CONTAMINATION_PREPARE.out.metrics,
        CONTAMINATION_SELECT_NUCLEOTIDE.out, nucleotide, ntDb, bam, taxdump,
        params.taxid, params.host_rank, manifest, helpers)
    CONTAMINATION_DIAGNOSTICS(CONTAMINATION_REPORT.out, plotter)
    CONTAMINATION_BLOBDIR(fasta, CONTAMINATION_REPORT.out, !!params.coverage_bam, taxdump, params.taxid, helpers)
}
