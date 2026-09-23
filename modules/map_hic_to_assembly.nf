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
    tuple val(meta), path(assembly_fasta), path(hic_r1, stageAs: 'hic_r1??/*'), path(hic_r2, stageAs: 'hic_r2??/*'), val(stage)

    output:
    tuple val(meta), val(stage),
          path("${meta.id}.sorted.bam"),
          //path("${meta.id}.sorted.bam.bai"),
          path("${meta.id}.sorted.bam.csi"),
          emit: bam

    tuple val(meta), val(stage),
          path("${meta.id}_mapping_stats.txt"),
          emit: stats

    path "versions.tsv", emit: versions
    tuple val(meta), val(stage), path("${meta.id}.readsets.tsv"), emit: readsets
    tuple val(meta), val(stage), path("${meta.id}.readset*.mapping_stats.txt"), emit: readset_stats

    script:
    def extra_args = (params.bwa_mem2_hic_args ?: "").toString()
    // samtools sort memory, scaled to the SLURM allocation:
    // ~1/4 of task.memory for sort buffers — bwa-mem2 (index + mapping, running concurrently
    // in the same pipe) gets the bulk — split across a capped sort-thread count.
    def sort_threads = Math.min(task.cpus as int, 16)
    def sort_mem_mb  = Math.max(768L, (task.memory.toMega() / 4 / sort_threads) as long)
    def r1s = hic_r1 instanceof Collection ? hic_r1 : [hic_r1]
    def r2s = hic_r2 instanceof Collection ? hic_r2 : [hic_r2]
    def sets = meta.hic_readsets
    if (!sets || sets.size() != r1s.size() || r1s.size() != r2s.size())
        error "${meta.id}: mapping read-set identities and mate lists do not agree"
    def bwa_threads = Math.max(1, (task.cpus as int) - 2)
    def mapping = sets.withIndex().collect { rs, i ->
        def prefix = String.format('readset%04d', i + 1)
        // RG identifies the LIBRARY, not the run: dedup may match across its runs.
        def rg = "@RG\\tID:${rs.library_id}\\tSM:${meta.sample}\\tLB:${rs.library_id}\\tPL:ILLUMINA"
        """
        bwa-mem2 mem -t ${bwa_threads} -5SP ${extra_args} -R '${rg}' '${assembly_fasta}' '${r1s[i]}' '${r2s[i]}' \\
          | awk 'BEGIN { FS=OFS="\\t" } /^@/ { print; next } { \$1="${prefix}_" \$1; if(length(\$1)>254) exit 2; print }' \\
          | samtools view -@ 1 -b -o ${prefix}.unsorted.bam -
        samtools sort -@ ${sort_threads} -m ${sort_mem_mb}M -T "\$PWD/${prefix}.sorttmp" \\
          -o ${prefix}.sorted.bam ${prefix}.unsorted.bam
        rm -f ${prefix}.unsorted.bam
        samtools flagstat ${prefix}.sorted.bam > ${meta.id}.${prefix}.mapping_stats.txt
        """
    }.join('\n')
    def bamNames = (0..<sets.size()).collect { String.format('readset%04d.sorted.bam', it + 1) }
    def merge = sets.size() == 1 ? "mv ${bamNames[0]} ${meta.id}.sorted.bam"
        : "samtools merge -c -@ ${sort_threads} -o ${meta.id}.sorted.bam ${bamNames.join(' ')}"
    def manifest = 'sample_id\tlibrary_id\treadset_id\tqname_prefix\n' + sets.withIndex().collect { rs, i ->
        "${meta.sample}\t${rs.library_id}\t${rs.readset_id}\t${String.format('readset%04d_', i + 1)}"
    }.join('\n')
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
    ${mapping}
    ${merge}
    rm -f ${bamNames.join(' ')}
    cat > ${meta.id}.readsets.tsv <<'READSET_MANIFEST'
${manifest}
READSET_MANIFEST

    # samtools index -@ ${task.cpus} ${meta.id}.sorted.bam
    samtools index -@ ${task.cpus} -c ${meta.id}.sorted.bam

    # -------------------------------------------------------------------------
    # 2) Mapping QC
    # -------------------------------------------------------------------------
    samtools flagstat -@ ${task.cpus} \\
      ${meta.id}.sorted.bam \\
      > ${meta.id}_mapping_stats.txt

    printf 'bwa-mem2\t%s\n' "\$(bwa-mem2 version 2>&1 | grep -m1 -oE '[0-9]+[.][0-9][0-9.]*' || echo unknown)" > versions.tsv
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    #touch ${meta.id}.sorted.bam.bai
    touch ${meta.id}.sorted.bam.csi
    touch ${meta.id}_mapping_stats.txt
    touch versions.tsv
    touch ${meta.id}.readsets.tsv
    touch ${meta.id}.readset0001.mapping_stats.txt
    """
}
