/*
========================================================================================
    BUILD MERYL DATABASE MODULE
========================================================================================
    Repo location: modules/build_meryl_db.nf

    Builds a k-mer database for Merqury from the sample's QC reads. Run once per sample,
    reused across all assembly-QC steps. Carries the sample-level meta (meta.id == sample).

    Phase 4a-iii: read-source-aware. `reads` is whatever the caller selected for QC —
    a single HiFi FASTQ for HiFi samples, or the Illumina R1+R2 pair for short-read
    samples. `meryl count` accepts multiple input files, so a PE pair builds one combined
    read k-mer DB (the standard Merqury approach for short reads).
========================================================================================
*/

process BUILD_MERYL_DB {
    tag "${meta.id}"
    label 'merqury'

    //publishDir "${params.outdir}/meryl_db", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.meryl"), emit: meryl_db

    script:
    def k = params.merqury_k ?: 21
    """
    # Build k-mer database from the sample's QC reads (1 HiFi FASTQ, or R1+R2 for PE).
    meryl count \\
        k=${k} \\
        threads=${task.cpus} \\
        memory=${task.memory.toGiga()} \\
        ${reads} \\
        output ${meta.id}.meryl
    """

    stub:
    """
    mkdir -p ${meta.id}.meryl
    touch ${meta.id}.meryl/merylIndex
    """
}
