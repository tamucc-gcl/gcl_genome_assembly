nextflow.enable.dsl = 2

include { CONTAMINATION_VERIFY_PROTEIN; CONTAMINATION_PREPARE; CONTAMINATION_BLASTX; CONTAMINATION_SELECT_NUCLEOTIDE;
          CONTAMINATION_BLASTN as CONTAMINATION_MEGABLAST;
          CONTAMINATION_BLASTN as CONTAMINATION_TARGETED_BLASTN;
          CONTAMINATION_REPORT; CONTAMINATION_DIAGNOSTICS; CONTAMINATION_BLOBDIR
} from '../modules/contamination_evidence.nf'

// Include one empty result per assembly so all-ambiguous sequences and empty shortlists survive the join.
def hitsByAssembly(sources, hits) {
    sources.map { id, fasta, host -> tuple(id, []) }
        .mix(hits.map { id, hit -> tuple(id, [hit]) })
        .groupTuple()
        .map { id, batches -> tuple(id, batches.flatten().sort { it.name }) }
}

workflow CONTAMINATION_SEARCH {
    take:
    assemblies     // (id, fasta, host taxid)
    coverage       // (id, optional BAM, optional BAI, optional coverage TSV)
    taxdump        // value
    settings       // database paths, thresholds and publication scope
    base           // pipeline source directory
    nt_manifest    // value, optional managed database manifest

    main:
    def sources = assemblies
    def coverageInput = coverage
    def bundle = file(settings.protein, checkIfExists: true).toAbsolutePath()
    def database = file(bundle.resolve('reference_proteomes.dmnd'), checkIfExists: true).toString()
    def manifest = file(bundle.resolve('protein_database.json'), checkIfExists: true)
    def helpers = ['contamination_review.py', 'contamination_evidence.py'].collect {
        file(base.resolve("py_scripts/${it}"), checkIfExists: true)
    }
    def reader = file(base.resolve('py_scripts/contamination_read_support.py'), checkIfExists: true)
    def plotter = file(base.resolve('r_scripts/plot_contamination_diagnostics.R'), checkIfExists: true)
    def mega = settings.mega
    def targeted = mega && settings.targeted

    CONTAMINATION_VERIFY_PROTEIN(file(database), manifest, taxdump, helpers)
    CONTAMINATION_PREPARE(sources.map { id, fasta, host -> tuple(id, fasta) },
        settings.stage, settings.window, settings.chunk, helpers[0])
    def chunks = CONTAMINATION_PREPARE.out.prepared.flatMap { id, metrics, queries ->
        queries.resolve('chunks.tsv').readLines().drop(1).collect { line ->
            tuple(id, queries.resolve(line.split('\t', 2)[0]))
        }
    }
    CONTAMINATION_BLASTX(chunks, database, CONTAMINATION_VERIFY_PROTEIN.out)
    def megaHits = Channel.empty()
    if (mega) {
        CONTAMINATION_MEGABLAST(chunks, settings.nt, 'megablast')
        megaHits = CONTAMINATION_MEGABLAST.out
    }

    def prepared = sources.join(CONTAMINATION_PREPARE.out.prepared, failOnDuplicate: true)
    def selectionInput = prepared.join(hitsByAssembly(sources, CONTAMINATION_BLASTX.out), failOnDuplicate: true)
        .join(hitsByAssembly(sources, megaHits), failOnDuplicate: true)
        .map { id, fasta, host, metrics, queries, hits, nt -> tuple(id, fasta, queries, host, hits, nt) }
    CONTAMINATION_SELECT_NUCLEOTIDE(selectionInput, mega, settings.followup, taxdump, settings.rank, helpers)
    def nucleotide = Channel.empty()
    if (targeted) {
        def targets = CONTAMINATION_SELECT_NUCLEOTIDE.out.selected.flatMap { id, selected ->
            selected.toFile().listFiles().findAll { it.name ==~ /target_.*\.fasta/ }
                .sort { it.name }.collect { tuple(id, it.toPath()) }
        }
        CONTAMINATION_TARGETED_BLASTN(targets, settings.nt, 'blastn')
        nucleotide = CONTAMINATION_TARGETED_BLASTN.out
    }
    def reportInput = prepared.join(CONTAMINATION_SELECT_NUCLEOTIDE.out.selected, failOnDuplicate: true)
        .join(hitsByAssembly(sources, nucleotide), failOnDuplicate: true)
        .join(coverageInput, failOnDuplicate: true, failOnMismatch: true)
        .map { id, fasta, host, metrics, queries, selected, hits, bam, bai, cov ->
            tuple(id, metrics, queries, selected, hits, bam, bai, cov, host)
        }
    CONTAMINATION_REPORT(reportInput, settings.nt, targeted, nt_manifest,
        taxdump, settings.rank, manifest, helpers, settings.nested, reader)
    CONTAMINATION_DIAGNOSTICS(CONTAMINATION_REPORT.out, plotter, settings.nested)
    def blobInput = sources.join(CONTAMINATION_REPORT.out, failOnDuplicate: true)
        .join(coverageInput, failOnDuplicate: true, failOnMismatch: true)
        .map { id, fasta, host, review, bam, bai, cov -> tuple(id, fasta, review, !!(bam || cov), host) }
    CONTAMINATION_BLOBDIR(blobInput, taxdump, helpers, settings.nested)

    emit:
    review = CONTAMINATION_REPORT.out
}
