/*
========================================================================================
    FINALIZE ASSEMBLY MODULE
========================================================================================
    Takes the final assembly (from gap filling or whichever stage is last)
    and produces a maximally compressed .fasta.gz for distribution.

    Uses bgzip (from htslib) for block-gzip compression which is:
    - Compatible with samtools/htslib tools (random access with .gzi index)
    - Maximally compressed at the highest compression level
    - Decompressible with standard gzip/gunzip

    Input:
    - tuple(haplotype_id, assembly_fasta)

    Output:
    - tuple(haplotype_id, compressed_fasta_gz)  — the .fasta.gz file
    - tuple(haplotype_id, gzi_index)            — the .fasta.gz.gzi index
    - tuple(haplotype_id, fai_index)            — the .fasta.gz.fai index
========================================================================================
*/

process FINALIZE_ASSEMBLY {
    tag "${haplotype_id}"
    label 'finalize_assembly'

    publishDir "${params.outdir}/assembly/final", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(assembly_fasta)

    output:
    tuple val(haplotype_id), path("${haplotype_id}.fasta.gz"),     emit: assembly
    tuple val(haplotype_id), path("${haplotype_id}.fasta.gz.gzi"), emit: gzi
    tuple val(haplotype_id), path("${haplotype_id}.fasta.gz.fai"), emit: fai

    script:
    """
    set -euo pipefail

    # bgzip at max compression level (-l 9), using multiple threads
    # --stdout so we control the output filename
    bgzip -c -l 9 --threads ${task.cpus} ${assembly_fasta} > ${haplotype_id}.fasta.gz

    # Create block-gzip index (.gzi) for random access
    bgzip --reindex ${haplotype_id}.fasta.gz

    # Create fasta index (.fai) — samtools can index bgzipped fastas directly
    samtools faidx ${haplotype_id}.fasta.gz
    """

    stub:
    """
    touch ${haplotype_id}.fasta.gz
    touch ${haplotype_id}.fasta.gz.gzi
    touch ${haplotype_id}.fasta.gz.fai
    """
}