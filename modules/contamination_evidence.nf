def quote(value) {
    "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process CONTAMINATION_PREPARE {
    label 'contamination_evidence'
    conda 'conda-forge::python=3.11'
    memory '8 GB'
    time '4h'

    input:
    path fasta
    val assembly
    val stage
    val window_bp
    val chunk_bp
    val database
    path manifest
    val taxdump
    path helpers

    output:
    path 'metrics', emit: metrics
    path 'queries', emit: queries
    path 'queries/protein_*.fasta', optional: true, emit: chunks

    script:
    """
    python3 contamination_evidence.py verify-database --db ${quote(database)} \
        --database ${quote(manifest)} --taxdump ${quote(taxdump)}
    python3 contamination_review.py prepare --fasta ${quote(fasta)} \
        --assembly ${quote(assembly)} --stage ${quote(stage)} --out metrics
    python3 contamination_review.py search-input --fasta ${quote(fasta)} \
        --assembly ${quote(assembly)} --stage ${quote(stage)} --out queries \
        --window-bp ${window_bp} --chunk-bp ${chunk_bp}
    """
}

process CONTAMINATION_BLASTX {
    tag "${chunk.name}"
    conda 'bioconda::diamond=2.1.14'
    cpus 8
    memory '24 GB'
    time '12h'

    input:
    path chunk
    val database

    output:
    path "${chunk.baseName}.hits.tsv"

    script:
    """
    diamond blastx --query ${quote(chunk)} --db ${quote(database)} \
        --out ${quote(chunk.baseName + '.hits.tsv')} --threads ${task.cpus} \
        --sensitive --evalue 1e-5 --frameshift 15 --top 20 --range-culling \
        --outfmt 6 qseqid sseqid staxids qstart qend bitscore evalue pident
    """
}

process CONTAMINATION_SELECT_NUCLEOTIDE {
    conda 'conda-forge::python=3.11'
    memory '8 GB'
    time '4h'

    input:
    path fasta
    path queries
    path hits
    val taxdump
    val host
    val rank
    path helpers

    output:
    path 'selected'

    script:
    def hitArgs = hits instanceof List ? hits.collect { quote(it) }.join(' ') : quote(hits)
    """
    python3 contamination_evidence.py select --fasta ${quote(fasta)} --prepared ${quote(queries)} \
        --taxdump ${quote(taxdump)} --taxid ${host} --host-rank ${quote(rank)} \
        --out selected --hits ${hitArgs}
    """
}

process CONTAMINATION_BLASTN {
    conda 'bioconda::blast=2.16.0'
    cpus 8
    memory '16 GB'
    time '24h'

    input:
    path selected
    val nt_db

    output:
    path 'nucleotide.hits.tsv'

    script:
    """
    blastdbcmd -db ${quote(nt_db)} -info
    if [ -s ${quote(selected.resolve('targets.fasta'))} ]; then
        blastn -task blastn -db ${quote(nt_db)} -query ${quote(selected.resolve('targets.fasta'))} \
            -out nucleotide.hits.tsv -num_threads ${task.cpus} -evalue 1e-10 -max_target_seqs 25 \
            -outfmt '6 qseqid sseqid staxids qstart qend bitscore evalue pident'
    else
        touch nucleotide.hits.tsv
    fi
    """
}

process CONTAMINATION_REPORT {
    conda 'conda-forge::python=3.11 bioconda::samtools=1.21'
    memory '8 GB'
    time '8h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path queries
    path metrics
    path selected
    path nt_hits
    val nt_db
    path bam
    val taxdump
    val host
    val rank
    path manifest
    path helpers

    output:
    path 'review'

    script:
    def depth = bam ? "samtools coverage -q 20 --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY ${quote(bam)} > coverage.tsv" : ''
    def cov = bam ? '--coverage coverage.tsv' : ''
    def nt = nt_hits ? "--nucleotide-hits ${quote(nt_hits)} --nucleotide-db ${quote(nt_db)}" : ''
    """
    ${depth}
    python3 contamination_evidence.py summarize --prepared ${quote(queries)} \
        --metrics ${quote(metrics)} --selected ${quote(selected)} --database ${quote(manifest)} \
        --taxdump ${quote(taxdump)} --taxid ${host} --host-rank ${quote(rank)} \
        --out review ${cov} ${nt}
    """
}

process CONTAMINATION_DIAGNOSTICS {
    conda 'conda-forge::r-base=4.3.3 conda-forge::r-ggplot2=3.5.2'
    memory '4 GB'
    time '2h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path review
    path plotter

    output:
    path 'diagnostics'

    script:
    """
    Rscript ${quote(plotter)} ${quote(review.resolve('sequence_evidence.tsv'))} diagnostics
    """
}

process CONTAMINATION_BLOBDIR {
    conda 'bioconda::blobtoolkit=4.5.5'
    cpus 4
    memory '16 GB'
    time '8h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fasta
    path review
    val hasCoverage
    val taxdump
    val host
    path helpers

    output:
    path 'blobdir'

    script:
    def columns = 'sequence=identifiers,evidence=contamination_evidence' +
        (hasCoverage ? ',mean_depth=mean_read_depth' : '')
    def axes = hasCoverage ? '--key plot.y=mean_read_depth' : ''
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
