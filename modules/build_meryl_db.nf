/*
========================================================================================
    BUILD MERYL DATABASE MODULE
========================================================================================
    Repo location: modules/build_meryl_db.nf

    Builds a k-mer database from HiFi reads for Merqury. Run once per sample, reused
    across all assembly-QC steps. Carries the sample-level meta (meta.id == sample).
========================================================================================
*/

process BUILD_MERYL_DB {
    tag "${meta.id}"
    label 'merqury'

    //publishDir "${params.outdir}/meryl_db", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hifi_fastq)

    output:
    tuple val(meta), path("${meta.id}.meryl"), emit: meryl_db

    script:
    def k = params.merqury_k ?: 21
    """
    # Build k-mer database from HiFi reads
    meryl count \\
        k=${k} \\
        threads=${task.cpus} \\
        memory=${task.memory.toGiga()} \\
        ${hifi_fastq} \\
        output ${meta.id}.meryl
    """

    stub:
    """
    mkdir -p ${meta.id}.meryl
    touch ${meta.id}.meryl/merylIndex
    """
}
