def quote(value) {
    "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process CONTAMINATION_STAGED_PREPARE {
    tag "${assembly}"
    conda 'conda-forge::python=3.11 bioconda::samtools=1.21'
    memory '8 GB'
    time '8h'
    input:
    tuple val(assembly), path(fasta), path(bam), path(bai), path(coverage)
    val sampling
    path helpers
    output:
    tuple val(assembly), path('metrics'), path('queries'), emit: prepared
    tuple val(assembly), path(bam), path(bai), path('staged_coverage.tsv'), emit: coverage
    script:
    def depth = bam ? "samtools coverage -q 20 --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY ${quote(bam)} > staged_coverage.tsv" :
        (coverage ? "cp ${quote(coverage)} staged_coverage.tsv" : 'touch staged_coverage.tsv')
    """
    ${depth}
    python3 contamination_review.py prepare --fasta ${quote(fasta)} --assembly ${quote(assembly)} --stage post_fcs --out metrics
    python3 contamination_staged.py sample --fasta ${quote(fasta)} --metrics metrics --coverage staged_coverage.tsv \
        --out queries --window-bp ${sampling.window} --spacing-bp ${sampling.spacing} \
        --short-contig-bp ${sampling.shortContig} --gc-delta ${sampling.gcDelta} --depth-fold ${sampling.depthFold}
    """
}

process CONTAMINATION_STAGED_NUCLEOTIDE {
    tag "${assembly}"
    conda 'conda-forge::python=3.11'
    memory '8 GB'
    time '4h'
    input:
    tuple val(assembly), path(fasta), path(queries), val(host), path(hits)
    val taxdump
    val rank
    path helpers
    output:
    tuple val(assembly), path('nt_queries')
    script:
    def hitArgs = hits instanceof List ? hits.collect { quote(it) }.join(' ') : quote(hits)
    """
    python3 contamination_staged.py route --fasta ${quote(fasta)} --prepared ${quote(queries)} \
        --taxdump ${quote(taxdump)} --taxid ${host} --host-rank ${quote(rank)} \
        --out nt_queries --hits ${hitArgs}
    """
}

process CONTAMINATION_SELECT_STAGED {
    tag "${assembly}"
    conda 'conda-forge::python=3.11'
    memory '8 GB'
    time '4h'
    input:
    tuple val(assembly), path(fasta), path(queries), val(host), path(hits), path(mega_hits), path(nt_queries)
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
    def megaFlag = mega_enabled ? "--megablast --megablast-plan ${quote(nt_queries.resolve('windows.tsv'))}" : ''
    """
    python3 contamination_evidence.py select --fasta ${quote(fasta)} --prepared ${quote(queries)} \
        --taxdump ${quote(taxdump)} --taxid ${host} --host-rank ${quote(rank)} \
        --out selected --followup-max-bp ${followup_max_bp} ${megaFlag} \
        --hits ${hitArgs} --megablast-hits ${megaArgs}
    """
}
