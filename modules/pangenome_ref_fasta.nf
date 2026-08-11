/*
========================================================================================
    PANGENOME REFERENCE FASTA MODULE
========================================================================================
    Repo location: modules/pangenome_ref_fasta.nf

    Extract the reference sample's path sequence from the graph as a FASTA (+ .fai). The
    headers are the graph's PanSN reference-path names (<ref_name>#0#chrN_p), so downstream
    tools that surject graph alignments to BAM or run reference-based callers get a FASTA
    whose contig names match the graph. This is a graph product; read mapping/indexing is a
    separate (downstream) pipeline.

    Input : tuple(taxid, ref_name, gbz)
    Output: ref_fasta / ref_fai
========================================================================================
*/

process PANGENOME_REF_FASTA {
    tag "${taxid}"
    label 'cactus_tools'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(ref_name), path(gbz)

    output:
    tuple val(taxid), path("${taxid}.reference.fa"),     emit: ref_fasta
    tuple val(taxid), path("${taxid}.reference.fa.fai"), emit: ref_fai

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"

    vg paths -x ${gbz} -F -S ${ref_name} > ${taxid}.reference.fa
    samtools faidx ${taxid}.reference.fa
    """

    stub:
    """
    : > ${taxid}.reference.fa
    : > ${taxid}.reference.fa.fai
    """
}
