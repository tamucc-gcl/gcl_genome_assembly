def quote(value) {
    "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process CONTAMINATION_VERIFY_PROTEIN {
    conda 'conda-forge::python=3.11'
    memory '1 GB'
    time '2h'

    input:
    path database
    path manifest
    val taxdump
    path helpers

    output:
    path 'protein_database.verified'

    script:
    """
    python3 contamination_evidence.py verify-database --db ${quote(database)} \
        --database ${quote(manifest)} --taxdump ${quote(taxdump)}
    touch protein_database.verified
    """
}

process CONTAMINATION_PREPARE {
    tag "${assembly}"
    label 'contamination_evidence'
    conda 'conda-forge::python=3.11'
    memory '8 GB'
    time '4h'

    input:
    tuple val(assembly), path(fasta)
    val stage
    val window_bp
    val chunk_bp
    path helper

    output:
    tuple val(assembly), path('metrics'), path('queries'), emit: prepared

    script:
    """
    python3 contamination_review.py prepare --fasta ${quote(fasta)} \
        --assembly ${quote(assembly)} --stage ${quote(stage)} --out metrics
    python3 contamination_review.py search-input --fasta ${quote(fasta)} \
        --assembly ${quote(assembly)} --stage ${quote(stage)} --out queries \
        --window-bp ${window_bp} --chunk-bp ${chunk_bp}
    """
}

process CONTAMINATION_BLASTX {
    tag "${assembly}/${chunk.name}"
    maxForks (params.contamination_blastx_max_forks as int)
    conda 'bioconda::diamond=2.1.14'
    cpus 8
    memory '24 GB'
    time '12h'

    input:
    tuple val(assembly), path(chunk)
    val database
    path verification

    output:
    tuple val(assembly), path("${chunk.baseName}.hits.tsv")

    script:
    """
    diamond blastx --query ${quote(chunk)} --db ${quote(database)} \
        --out ${quote(chunk.baseName + '.hits.tsv')} --threads ${task.cpus} \
        --sensitive --evalue 1e-5 --frameshift 15 --top 20 --range-culling \
        --outfmt 6 qseqid sseqid staxids qstart qend bitscore evalue pident stitle
    """
}

process CONTAMINATION_SELECT_NUCLEOTIDE {
    tag "${assembly}"
    conda 'conda-forge::python=3.11'
    memory '8 GB'
    time '4h'

    input:
    tuple val(assembly), path(fasta), path(queries), val(host), path(hits), path(mega_hits)
    val mega_enabled
    val followup_max_bp
    val taxdump
    val rank
    path helpers

    output:
    tuple val(assembly), path('selected'), emit: selected

    script:
    def hitArgs = hits instanceof List ? hits.collect { quote(it) }.join(' ') : quote(hits)
    def megaArgs = mega_hits instanceof List ? mega_hits.collect { quote(it) }.join(' ') : quote(mega_hits)
    def megaFlag = mega_enabled ? '--megablast' : ''
    """
    python3 contamination_evidence.py select --fasta ${quote(fasta)} --prepared ${quote(queries)} \
        --taxdump ${quote(taxdump)} --taxid ${host} --host-rank ${quote(rank)} \
        --out selected --followup-max-bp ${followup_max_bp} ${megaFlag} \
        --hits ${hitArgs} --megablast-hits ${megaArgs}
    """
}

process CONTAMINATION_BLASTN {
    tag "${assembly}/${mode}/${chunk.name}"
    maxForks (params.contamination_blastn_max_forks as int)
    conda 'bioconda::blast=2.16.0'
    cpus 8
    memory '16 GB'
    time '24h'

    input:
    tuple val(assembly), path(chunk)
    val nt_db
    val mode

    output:
    tuple val(assembly), path("${mode}_${chunk.baseName}.hits.tsv")

    script:
    """
    blastn -task ${mode} -db ${quote(nt_db)} -query ${quote(chunk)} \
        -out ${quote(mode + '_' + chunk.baseName + '.hits.tsv')} -num_threads ${task.cpus} \
        -evalue 1e-10 -culling_limit 10 \
        -outfmt '6 qseqid saccver staxids qstart qend bitscore evalue pident stitle'
    """
}

process CONTAMINATION_REPORT {
    tag "${assembly}"
    conda 'conda-forge::python=3.11 bioconda::samtools=1.21 bioconda::pysam=0.23.3'
    memory '8 GB'
    time '8h'
    publishDir { nested ? "${params.outdir}/assembly/contamination/${assembly}" : params.outdir }, mode: 'copy'

    input:
    tuple val(assembly), path(metrics), path(queries), path(selected), path(nt_hits), path(bam), path(bai),
          path(coverage_table), val(host)
    val nt_db
    val targeted
    path nt_manifest
    val taxdump
    val rank
    path manifest
    path helpers
    val nested
    path reader

    output:
    tuple val(assembly), path('review')

    script:
    def depth = bam ? "samtools coverage -q 20 --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY ${quote(bam)} > coverage.tsv" : ''
    def cov = bam ? '--coverage coverage.tsv' :
        (coverage_table ? "--coverage ${quote(coverage_table)} --coverage-source provided-coverage-table" : '')
    def ntArgs = nt_hits instanceof List ? nt_hits.collect { quote(it) }.join(' ') : quote(nt_hits)
    def nt = (nt_db ? "--nucleotide-db ${quote(nt_db)}" : '') +
        (nt_manifest ? " --nucleotide-manifest ${quote(nt_manifest)}" : '')
    def targetedFlag = targeted ? '--targeted' : ''
    def readSupport = bai ? "python3 ${quote(reader)} --bam ${quote(bam)} --index ${quote(bai)} --review review" : ''
    """
    ${depth}
    python3 contamination_evidence.py summarize --prepared ${quote(queries)} \
        --metrics ${quote(metrics)} --selected ${quote(selected)} --database ${quote(manifest)} \
        --taxdump ${quote(taxdump)} --taxid ${host} --host-rank ${quote(rank)} \
        --out review ${cov} ${nt} ${targetedFlag} --nucleotide-hits ${ntArgs}
    ${readSupport}
    """
}

process CONTAMINATION_REBUILD_REPORT {
    conda 'conda-forge::python=3.11'
    memory '8 GB'
    time '8h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path previous, stageAs: 'prior'
    path fasta
    path helpers
    val assembly

    output:
    tuple val(assembly), path('review')

    script:
    """
    python3 contamination_evidence.py rebuild-report --previous ${quote(previous)} \
        --fasta ${quote(fasta)} --out review
    """
}

process CONTAMINATION_DIAGNOSTICS {
    conda 'conda-forge::r-base=4.3.3 conda-forge::r-ggplot2=3.5.2 conda-forge::r-patchwork=1.3.0'
    memory '4 GB'
    time '2h'
    publishDir { nested ? "${params.outdir}/assembly/contamination/${assembly}" : params.outdir }, mode: 'copy'

    input:
    tuple val(assembly), path(review)
    path plotter
    val nested

    output:
    path 'diagnostics'

    script:
    """
    Rscript ${quote(plotter)} ${quote(review.resolve('sequence_evidence.tsv'))} diagnostics \
        ${quote(review.resolve('locus_plot_hits.tsv'))} ${quote(review)}
    """
}

process CONTAMINATION_BLOBDIR {
    conda 'bioconda::blobtoolkit=4.5.5'
    cpus 4
    memory '16 GB'
    time '8h'
    publishDir { nested ? "${params.outdir}/assembly/contamination/${assembly}" : params.outdir }, mode: 'copy'

    input:
    tuple val(assembly), path(fasta), path(review), val(hasCoverage), val(host)
    val taxdump
    path helpers
    val nested

    output:
    path 'blobdir'

    script:
    def columns = 'sequence=identifiers,evidence=contamination_evidence,plot_taxon=screen_taxon' +
        (hasCoverage ? ',mean_depth=mean_read_depth' : '')
    def axes = '--key plot.cat=screen_taxon' + (hasCoverage ? ' --key plot.y=mean_read_depth' : '')
    """
    python3 contamination_evidence.py blob-input --fasta ${quote(fasta)} \
        --exclusions ${quote(review.resolve('blobdir_exclusions.tsv'))} --out blob-input.fasta
    if [ ! -s blob-input.fasta ]; then
        mkdir blobdir
        echo 'BlobToolKit requires at least one sequence with measurable GC. See ../review/blobdir_exclusions.tsv.' > blobdir/README.md
        exit 0
    fi
    python3 contamination_evidence.py blob-taxonomy --taxdump ${quote(taxdump)} --taxid ${host} \
        --hits ${quote(review.resolve('blob_hits.tsv'))} --out taxonomy
    blobtools create --fasta blob-input.fasta --taxid ${host} --taxdump taxonomy blobdir
    if [ -s taxonomy/hits.tsv ]; then
        blobtools add --hits taxonomy/hits.tsv \
            --hits-cols 1=qseqid,2=staxids,3=bitscore,4=sseqid,5=evalue,6=qstart,7=qend \
            --taxrule bestsum --taxdump taxonomy blobdir
    fi
    blobtools add --text ${quote(review.resolve('blob_sequence_evidence.tsv'))} \
        --text-delimiter '\\t' --text-header --text-cols ${quote(columns)} ${axes} blobdir
    """
}
