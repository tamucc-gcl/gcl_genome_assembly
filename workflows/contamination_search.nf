nextflow.enable.dsl = 2

include { CONTAMINATION_VERIFY_PROTEIN; CONTAMINATION_PREPARE; CONTAMINATION_BLASTX; CONTAMINATION_SELECT_NUCLEOTIDE;
          CONTAMINATION_BLASTN as CONTAMINATION_MEGABLAST;
          CONTAMINATION_BLASTN as CONTAMINATION_TARGETED_BLASTN;
          CONTAMINATION_REPORT; CONTAMINATION_DIAGNOSTICS; CONTAMINATION_BLOBDIR
} from '../modules/contamination_evidence.nf'
include { CONTAMINATION_MEGABLAST_NODES } from './contamination_megablast_nodes.nf'

// GString metadata and plain String IDs from JSON/cache must use identical join keys.
def assemblyKeys(channel) {
    channel.map { row -> [row[0].toString()] + row.drop(1) }
}

// Include one empty result per assembly so all-ambiguous sequences and empty shortlists survive the join.
def hitsByAssembly(counts, hits) {
    def keyedCounts = assemblyKeys(counts)
    keyedCounts.map { id, count -> tuple(groupKey(id, count + 1), []) }
        .mix(assemblyKeys(hits).combine(keyedCounts, by: 0)
            .map { id, hit, count -> tuple(groupKey(id, count + 1), [hit]) })
        .groupTuple()
        .map { key, batches -> tuple(key.getGroupTarget(), batches.flatten().sort { it.name }) }
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
    def sources = assemblyKeys(assemblies)
    def coverageInput = assemblyKeys(coverage)
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
    if (!(params.contamination_megablast_execution in ['chunks', 'nodes']))
        error 'Set --contamination_megablast_execution to chunks or nodes'
    ['contamination_megablast_nodes', 'contamination_megablast_slots', 'contamination_megablast_threads'].each { name ->
        if (!(params[name]?.toString() ==~ /[1-9][0-9]*/))
            error "Supply a positive integer for --${name}"
    }

    CONTAMINATION_VERIFY_PROTEIN(file(database), manifest, taxdump, helpers)
    CONTAMINATION_PREPARE(sources.map { id, fasta, host -> tuple(id, fasta) },
        settings.stage, settings.window, settings.chunk, helpers[0])
    def preparedOutputs = assemblyKeys(CONTAMINATION_PREPARE.out.prepared)
    def chunks = preparedOutputs.flatMap { id, metrics, queries ->
        queries.resolve('chunks.tsv').readLines().drop(1).collect { line ->
            tuple(id, queries.resolve(line.split('\t', 2)[0]))
        }
    }
    def counts = preparedOutputs.map { id, metrics, queries ->
        tuple(id, queries.resolve('chunks.tsv').readLines().size() - 1)
    }
    CONTAMINATION_BLASTX(chunks, database, CONTAMINATION_VERIFY_PROTEIN.out)
    def megaHits = Channel.empty()
    if (mega) {
        if (params.contamination_megablast_execution == 'nodes') {
            CONTAMINATION_MEGABLAST_NODES(preparedOutputs, settings.nt, base)
            megaHits = CONTAMINATION_MEGABLAST_NODES.out.results
        } else {
            CONTAMINATION_MEGABLAST(chunks, settings.nt, 'megablast')
            megaHits = CONTAMINATION_MEGABLAST.out
        }
    }

    def prepared = sources.join(preparedOutputs, failOnDuplicate: true, failOnMismatch: true)
    def selectionInput = prepared.join(hitsByAssembly(counts, CONTAMINATION_BLASTX.out), failOnDuplicate: true)
        .join(hitsByAssembly(mega ? counts : counts.map { id, count -> tuple(id, 0) }, megaHits), failOnDuplicate: true)
        .map { id, fasta, host, metrics, queries, hits, nt -> tuple(id, fasta, queries, host, hits, nt) }
    CONTAMINATION_SELECT_NUCLEOTIDE(selectionInput, mega, settings.followup, taxdump, settings.rank, helpers)
    def selectedOutputs = assemblyKeys(CONTAMINATION_SELECT_NUCLEOTIDE.out.selected)
    def nucleotide = Channel.empty()
    if (targeted) {
        def targets = selectedOutputs.flatMap { id, selected ->
            selected.toFile().listFiles().findAll { it.name ==~ /target_.*\.fasta/ }
                .sort { it.name }.collect { tuple(id, it.toPath()) }
        }
        CONTAMINATION_TARGETED_BLASTN(targets, settings.nt, 'blastn')
        nucleotide = CONTAMINATION_TARGETED_BLASTN.out
    }
    def targetCounts = selectedOutputs.map { id, selected ->
        tuple(id, targeted ? selected.toFile().listFiles().count { it.name ==~ /target_.*\.fasta/ } : 0)
    }
    def reportInput = prepared.join(selectedOutputs, failOnDuplicate: true)
        .join(hitsByAssembly(targetCounts, nucleotide), failOnDuplicate: true)
        .join(coverageInput, failOnDuplicate: true, failOnMismatch: true)
        .map { id, fasta, host, metrics, queries, selected, hits, bam, bai, cov ->
            tuple(id, metrics, queries, selected, hits, bam, bai, cov, host)
        }
    CONTAMINATION_REPORT(reportInput, settings.nt, targeted, nt_manifest,
        taxdump, settings.rank, manifest, helpers, settings.nested, reader)
    def reports = assemblyKeys(CONTAMINATION_REPORT.out)
    CONTAMINATION_DIAGNOSTICS(reports, plotter, settings.nested)
    def blobInput = sources.join(reports, failOnDuplicate: true)
        .join(coverageInput, failOnDuplicate: true, failOnMismatch: true)
        .map { id, fasta, host, review, bam, bai, cov -> tuple(id, fasta, review, !!(bam || cov), host) }
    CONTAMINATION_BLOBDIR(blobInput, taxdump, helpers, settings.nested)

    emit:
    review = reports
}
