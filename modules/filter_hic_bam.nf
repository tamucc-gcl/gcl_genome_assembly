/*
========================================================================================
    FILTER HI-C BAM MODULE (FAST + EXPLICIT)
========================================================================================
    Repo location: modules/filter_hic_bam.nf

    Filters Hi-C BAM to keep valid UU pairs, MAPQ filter at parse, PCR-dup removal in
    pair-space. Produces sorted/indexed BAM, filtered pairs.gz, and a combined stats report.
    Work-dir-local temp (job-namespaced via meta.id) avoids /tmp collisions across tasks.
    Stage parameter controls publishDir but not filenames.
========================================================================================
*/

process FILTER_HIC_BAM {
    tag "${meta.id}_${stage}"
    label 'filter_hic_bam'

    publishDir "${params.outdir}/bam/hic/${stage}/filtered", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(stage), path(bam), path(bai), path(assembly_fasta)

    output:
    tuple val(meta), val(stage), path("${meta.id}.filtered.sorted.bam"), path("${meta.id}.filtered.sorted.bam.bai"), emit: bam
    tuple val(meta), val(stage), path("${meta.id}_filtering_stats.txt"), emit: stats
    tuple val(meta), val(stage), path("${meta.id}.pairs.gz"), emit: pairs
    tuple val(meta), val(stage), path("${meta.id}_parse_stats.txt"), emit: parse_stats
    tuple val(meta), val(stage), path("${meta.id}_dedup_stats.txt"), emit: dedup_stats

    script:
    """
    set -euo pipefail
    export LC_ALL=C

    MINQ=${params.hic_min_mapq ?: 30}

    # Force all temp file defaults onto the scratch-backed working directory.
    # Nextflow's scratch directive already places \$PWD on local NVMe; this
    # ensures pairtools sort (which reads TMPDIR) doesn't fall back to /tmp.
    export TMPDIR="\$PWD"

    # Thread budget — divide task.cpus across concurrent pipe stages to avoid
    # over-subscription.  The sort stages are the bottleneck, so they get the
    # largest share.
    CPUS=${task.cpus}
    T_SORT=\$(( CPUS / 2 ))           # pairtools sort / samtools sort
    T_IO=\$(( (CPUS - T_SORT) / 2 ))  # collate, parse, split, view
    T_IO=\$(( T_IO > 1 ? T_IO : 2 ))  # minimum 2

    # -------------------------------------------------------------------------
    # 1) Prepare chrom sizes (prefer existing .fai if present)
    # -------------------------------------------------------------------------
    if [[ ! -s "${assembly_fasta}.fai" ]]; then
      samtools faidx ${assembly_fasta}
    fi
    cut -f1,2 ${assembly_fasta}.fai > chrom.sizes

    # -------------------------------------------------------------------------
    # 2) Parse -> sort -> dedup -> select (keep only valid UU pairs)
    #    IMPORTANT: output is a .pairsam.gz containing sam1/sam2 so we can restore BAM.
    # -------------------------------------------------------------------------
    samtools collate -T "\$PWD/${meta.id}_collate" -@ \${T_IO} -O -u ${bam} | \\
      pairtools parse \\
        --min-mapq \${MINQ} \\
        --walks-policy 5unique \\
        --max-inter-align-gap 30 \\
        --chroms-path chrom.sizes \\
        --output-stats ${meta.id}_parse_stats.txt \\
        --nproc-in  \${T_IO} \\
        --nproc-out \${T_IO} \\
        - \\
    | pairtools sort \\
        --nproc \${T_SORT} \\
        --tmpdir "\$PWD" \\
    | pairtools dedup \\
        --mark-dups \\
        --output-stats ${meta.id}_dedup_stats.txt \\
        --output-dups ${meta.id}.dups.pairs.gz \\
    | pairtools select '(pair_type == "UU")' \\
        --output ${meta.id}.pairsam.gz

    # -------------------------------------------------------------------------
    # 3) Create filtered .pairs.gz AND restore a filtered BAM, fast:
    #    pairtools split streams SAM to stdout; samtools consumes + sorts once.
    # -------------------------------------------------------------------------
    pairtools split \\
      --output-pairs ${meta.id}.pairs.gz \\
      --output-sam - \\
      --nproc-in  \${T_IO} \\
      --nproc-out \${T_IO} \\
      ${meta.id}.pairsam.gz \\
    | samtools view -@ \${T_IO} -b - \\
    | samtools sort -@ \${T_SORT} -T "\$PWD/${meta.id}.sort" -o ${meta.id}.filtered.sorted.bam -

    samtools index -@ \${CPUS} ${meta.id}.filtered.sorted.bam

    # -------------------------------------------------------------------------
    # 4) Optional pair-level stats (pairtools stats) if available
    # -------------------------------------------------------------------------
    if pairtools stats --help >/dev/null 2>&1; then
      pairtools stats ${meta.id}.pairs.gz > ${meta.id}_pairs_stats.txt
      HAVE_PAIRTOOLS_STATS=1
    else
      HAVE_PAIRTOOLS_STATS=0
    fi

    # -------------------------------------------------------------------------
    # 5) Comprehensive statistics report
    #    NOTE: BAM counts are alignment-level; pairs.gz / pairtools stats are pair-level.
    # -------------------------------------------------------------------------
    # Exclude secondary (0x100) and supplementary (0x800) alignments:
    # 0x100 + 0x800 = 0x900 = 2304
    original_primary=\$(samtools view -c -F 2304 ${bam})
    filtered_primary=\$(samtools view -c -F 2304 ${meta.id}.filtered.sorted.bam)

    # Avoid hard dependency on bc; awk is everywhere
    retained_pct=\$(awk -v f="\$filtered_primary" -v o="\$original_primary" 'BEGIN{ if(o>0) printf("%.2f", (f/o)*100); else print "NA"; }')

    {
      echo "================================================================================"
      echo "Hi-C BAM Filtering Statistics for ${meta.id} (${stage})"
      echo "Generated: \$(date)"
      echo "pairtools version: \$(pairtools --version 2>/dev/null || echo 'unknown')"
      echo "samtools version:  \$(samtools --version | head -n 1)"
      echo "================================================================================"
      echo
      echo "Filtering definition (explicit):"
      echo "  - pairtools parse: min MAPQ = \${MINQ}"
      echo "  - pairtools dedup: PCR duplicates marked/removed in pair-space"
      echo "  - pairtools select: keep only UU pairs (both ends uniquely mapped)"
      echo
      echo "--------------------------------------------------------------------------------"
      echo "ORIGINAL BAM (samtools flagstat)"
      echo "--------------------------------------------------------------------------------"
      samtools flagstat ${bam}
      echo
      echo "--------------------------------------------------------------------------------"
      echo "PAIRTOOLS PARSE STATS"
      echo "--------------------------------------------------------------------------------"
      cat ${meta.id}_parse_stats.txt
      echo
      echo "--------------------------------------------------------------------------------"
      echo "PAIRTOOLS DEDUP STATS"
      echo "--------------------------------------------------------------------------------"
      cat ${meta.id}_dedup_stats.txt
      echo
      if [[ "\${HAVE_PAIRTOOLS_STATS}" -eq 1 ]]; then
        echo "--------------------------------------------------------------------------------"
        echo "PAIRTOOLS PAIRS STATS (filtered pairs.gz)"
        echo "--------------------------------------------------------------------------------"
        cat ${meta.id}_pairs_stats.txt
        echo
      fi
      echo "--------------------------------------------------------------------------------"
      echo "FILTERED BAM (samtools flagstat)"
      echo "--------------------------------------------------------------------------------"
      samtools flagstat ${meta.id}.filtered.sorted.bam
      echo
      echo "--------------------------------------------------------------------------------"
      echo "ALIGNMENT-LEVEL RETENTION (primary alignments only; excludes secondary+supplementary)"
      echo "--------------------------------------------------------------------------------"
      echo "Original primary alignments:  \${original_primary}"
      echo "Filtered primary alignments:  \${filtered_primary}"
      echo "Retention rate:              \${retained_pct}%"
      echo
      echo "Notes:"
      echo "  - Use pairs.gz / pairtools stats for true pair-level retention."
      echo "  - UU-only is typically a good default for scaffolding to reduce noise."
      echo "================================================================================"
    } > ${meta.id}_filtering_stats.txt

    # Cleanup intermediates
    rm -f ${meta.id}.pairsam.gz
    rm -f ${meta.id}.dups.pairs.gz
    rm -f chrom.sizes
    """

    stub:
    """
    touch ${meta.id}.filtered.sorted.bam
    touch ${meta.id}.filtered.sorted.bam.bai
    touch ${meta.id}_filtering_stats.txt
    touch ${meta.id}.pairs.gz
    touch ${meta.id}_parse_stats.txt
    touch ${meta.id}_dedup_stats.txt
    """
}
