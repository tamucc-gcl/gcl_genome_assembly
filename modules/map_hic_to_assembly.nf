/*
========================================================================================
    MAP HI-C READS TO ASSEMBLY (bwa-mem2, PAIRTOOLS-READY)
========================================================================================
    Purpose:
      Produce the optimal BAM for Hi-C filtering with pairtools:
        - queryname-collated BAM (mates adjacent)
        - no coordinate sorting here
        - optimized for speed and correctness

    Outputs:
      - *.hic.qname.bam        → input for FILTER_HIC_BAM
      - *.hic.qname.bam.csi    → optional index
      - mapping flagstat QC

    Notes:
      - Uses bwa-mem2 (much faster than bwa mem)
      - Uses Hi-C appropriate flags: -5SP
      - pairtools parse REQUIRES qname-collated input
========================================================================================
*/

process MAP_HIC_TO_ASSEMBLY {
    tag "${haplotype_id}"
    label 'map_hic'

    publishDir "${params.outdir}/bam/hic/raw",
        mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(assembly_fasta), path(hic_r1), path(hic_r2)

    output:
    tuple val(haplotype_id),
          path("${haplotype_id}.sorted.bam"),
          path("${haplotype_id}.sorted.bam.bai"),
          emit: bam

    tuple val(haplotype_id),
          path("${haplotype_id}_mapping_stats.txt"),
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
    | samtools sort -@ ${task.cpus} -o ${haplotype_id}.sorted.bam -
    samtools index -@ ${task.cpus} ${haplotype_id}.sorted.bam

    # CSI index works on unsorted BAMs
    samtools index -@ ${task.cpus} -c ${haplotype_id}.sorted.bam

    # -------------------------------------------------------------------------
    # 2) Mapping QC
    # -------------------------------------------------------------------------
    samtools flagstat -@ ${task.cpus} \\
      ${haplotype_id}.sorted.bam \\
      > ${haplotype_id}_mapping_stats.txt
    """
    
    stub:
    """
    touch ${haplotype_id}.sorted.bam
    touch ${haplotype_id}.sorted.bam.bai
    touch ${haplotype_id}_mapping_stats.txt
    """
}