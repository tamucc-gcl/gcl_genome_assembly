process MAP_HIFI_FINAL {
    tag "${meta.id}"
    label 'map_hifi_final'
    publishDir "${params.outdir}/bam/hifi/final", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly), path(hifi_reads)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.reference.fai"), emit: fai
    tuple val(meta), path("${meta.id}.reference.sha256"), emit: reference_digest
    path 'versions.tsv', emit: versions

    script:
    def map_threads = Math.max(1, (task.cpus as int) - 2)
    """
    set -euo pipefail
    samtools faidx ${assembly}
    cp ${assembly}.fai ${meta.id}.reference.fai
    sha256sum ${assembly} > ${meta.id}.reference.sha256
    minimap2 -t ${map_threads} -ax map-hifi ${assembly} ${hifi_reads} \\
      | samtools sort -@ 1 -m 1G -o ${meta.id}.sorted.bam
    samtools quickcheck ${meta.id}.sorted.bam
    samtools index -@ 1 ${meta.id}.sorted.bam
    printf 'minimap2\\t%s\\n' "\$(minimap2 --version)" > versions.tsv
    printf 'samtools\\t%s\\n' "\$(samtools --version | head -n1)" >> versions.tsv
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    touch ${meta.id}.reference.fai ${meta.id}.reference.sha256 versions.tsv
    """
}
