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

    path "versions.tsv", emit: versions

    script:
    def extra_args = (params.bwa_mem2_hic_args ?: "").toString()
    // samtools sort memory, scaled to the SLURM allocation:
    // ~1/4 of task.memory for sort buffers — bwa-mem2 (index + mapping, running concurrently
    // in the same pipe) gets the bulk — split across a capped sort-thread count.
    def sort_threads = Math.min(task.cpus as int, 16)
    def sort_mem_mb  = Math.max(768L, (task.memory.toMega() / 4 / sort_threads) as long)
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
    # Map to an unsorted BAM first, then sort as a separate step.
    # (A fused bwa|view|sort pipe keeps sort running under bwa's backpressure for the
    #  whole mapping run; on large Hi-C sets that long-lived pipe dies mid-merge -> SIGPIPE/141.)
    bwa-mem2 mem -t ${task.cpus} -5SP ${extra_args} ${assembly_fasta} ${hic_r1} ${hic_r2} \
      | samtools view -@ ${task.cpus} -b -o ${meta.id}.unsorted.bam -
    samtools sort -@ ${sort_threads} -m ${sort_mem_mb}M -T "\$PWD/${meta.id}.sorttmp" \
      -o ${meta.id}.sorted.bam ${meta.id}.unsorted.bam
    rm -f ${meta.id}.unsorted.bam

    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam

    # CSI index works on unsorted BAMs
    # samtools index -@ ${task.cpus} -c ${meta.id}.sorted.bam

    # -------------------------------------------------------------------------
    # 2) Mapping QC
    # -------------------------------------------------------------------------
    samtools flagstat -@ ${task.cpus} \\
      ${meta.id}.sorted.bam \\
      > ${meta.id}_mapping_stats.txt

    # bwa-mem2's launcher prints a multi-line banner before the version; a bare `... | head`
    # closes the pipe after line 1 and SIGPIPEs bwa-mem2 -> exit 141 under pipefail. Wrapping in
    # $() contains the SIGPIPE (as the other version captures do); the regex grabs the real
    # version token, not the banner (and skips the "512" in avx512bw by requiring a dot).
    printf 'bwa-mem2\t%s\n' "\$(bwa-mem2 version 2>&1 | grep -m1 -oE '[0-9]+[.][0-9][0-9.]*' || echo unknown)" > versions.tsv
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    touch ${meta.id}.sorted.bam.bai
    touch ${meta.id}_mapping_stats.txt
    touch versions.tsv
    """
}
