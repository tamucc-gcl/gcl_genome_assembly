/*
========================================================================================
    BUILD MERYL DATABASE MODULE
========================================================================================
    Builds k-mer database from HiFi reads for use with Merqury
    - Run once per sample
    - Reused across all assembly QC steps
    - Significantly speeds up iterative QC
========================================================================================
*/

process BUILD_MERYL_DB {
    tag "${sample_id}"
    label 'merqury'
    
    //publishDir "${params.outdir}/meryl_db", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hifi_fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}.meryl"), emit: meryl_db
    
    script:
    def k = params.merqury_k ?: 21
    """
    # Build k-mer database from HiFi reads
    meryl count \\
        k=${k} \\
        threads=${task.cpus} \\
        memory=${task.memory.toGiga()} \\
        ${hifi_fastq} \\
        output ${sample_id}.meryl
    """
    
    stub:
    """
    mkdir -p ${sample_id}.meryl
    touch ${sample_id}.meryl/merylIndex
    """
}