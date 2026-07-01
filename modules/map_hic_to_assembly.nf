/*
========================================================================================
    MAP HI-C READS TO ASSEMBLY (bwa-mem2, PAIRTOOLS-READY)
========================================================================================
    Repo location: modules/map_hic_to_assembly.nf

    Produces a sorted, indexed BAM for Hi-C filtering with pairtools.
    Notes:
      - bwa-mem2, Hi-C flags -5SP
      - Stage parameter controls publishDir but not filenames
========================================================================================
*/

process MAP_HIC_TO_ASSEMBLY {
    tag "${meta.id}_${stage}"
    label 'map_hic'

    publishDir "${params.outdir}/bam/hic/${stage}/raw",
        mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly_fasta), path(hic_r1), path(hic_r2), val(stage)

    output:
    tuple val(meta), val(stage),
          path("${meta.id}.sorted.bam"),
          path("${meta.id}.sorted.bam.bai"),
          emit: bam

    tuple val(meta), val(stage),
          path("${meta.id}_mapping_stats.txt"),
          emit: stats

    script:
    def extra_args = (params.bwa_mem2_hic_args ?: "").toString()
    """
    set -euo pipefail
    export LC_ALL=C
    TMPDIR="\${TMPDIR:-\$PWD}"

    # -------------------------------------------------------------------------
    # 0) Reference indexing (guarded)
    #    bwa-mem2 index creates .0123/.amb/.ann/.bwt.2bit
    # -------------------------------------------------------------------------
    if [[ ! -s "${assembly_fasta}.fai" ]]; then
      samtools faidx ${assembly_fasta}
    fi

    if [[ ! -s "${assembly_fasta}.0123" ]]; then
      bwa-mem2 index ${assembly_fasta}
    fi

    # -------------------------------------------------------------------------
    # 1) Map Hi-C reads + collate by queryname
    #    samtools collate keeps mates together for pairtools
    # -------------------------------------------------------------------------
    bwa-mem2 mem -t ${task.cpus} -5SP ${extra_args} ${assembly_fasta} ${hic_r1} ${hic_r2} \
    | samtools view -@ ${task.cpus} -b - \
    | samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam -
    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam

    # CSI index works on unsorted BAMs
    samtools index -@ ${task.cpus} -c ${meta.id}.sorted.bam

    # -------------------------------------------------------------------------
    # 2) Mapping QC
    # -------------------------------------------------------------------------
    samtools flagstat -@ ${task.cpus} \\
      ${meta.id}.sorted.bam \\
      > ${meta.id}_mapping_stats.txt
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    touch ${meta.id}.sorted.bam.bai
    touch ${meta.id}_mapping_stats.txt
    """
}
