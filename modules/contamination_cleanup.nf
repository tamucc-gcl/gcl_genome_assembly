process CONTAMINATION_DECISIONS {
    conda 'conda-forge::python=3.11'
    input:
    path sheet
    val assemblies
    path helper

    output:
    path 'decisions.json'

    script:
    def ids = assemblies.collect { "'"+it.toString().replace("'", "'\"'\"'")+"'" }.join(' ')
    """
    python3 contamination_review.py plan --review "${sheet}" --out decisions.json --assemblies ${ids}
    """
}

process APPLY_CONTAMINATION_REVIEW {
    tag "${meta.id}"
    conda 'conda-forge::python=3.11'
    memory '6 GB'
    time '4h'
    publishDir "${params.outdir}/assembly/contig/decontam/blobtools/${meta.id}/cleanup", mode: 'copy', pattern: 'audit'

    input:
    tuple val(meta), path(fasta)
    path sheet
    path helper

    output:
    tuple val(meta), path('audit/retained.fasta'), emit: assembly
    path 'audit'

    script:
    """
    python3 contamination_review.py apply --fasta "${fasta}" --assembly "${meta.id}" \
        --stage post_fcs --review "${sheet}" --out audit
    """
}
