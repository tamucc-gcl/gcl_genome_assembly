/*
========================================================================================
    FILTER HI-C BAM MODULE (FAST + EXPLICIT)
========================================================================================
    Filters Hi-C BAM files to remove invalid/low-quality pairs and PCR duplicates
    - Pair-level filtering (explicit): keeps only UU pairs (both ends uniquely mapped)
    - MAPQ filter at parse stage (default: 30)
    - PCR duplicate removal using pairtools (in pair-space)
    - Produces:
        1) coordinate-sorted, indexed BAM (for scaffolding / visualization)
        2) filtered pairs.gz (for QC/contact-map tooling)
        3) stats report combining samtools + pairtools summaries

    Speed improvements vs original:
    - Avoids writing gigantic intermediate SAMs and merge steps
    - Streams pairsam -> samtools sort directly
    - Uses TMPDIR for sorting scratch space
    - Uses pairtools multithreading where supported
========================================================================================
*/

process FILTER_HIC_BAM {
    tag "${haplotype_id}"
    label 'filter_hic_bam'

    publishDir "${params.outdir}/bam/hic/filtered", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(bam), path(bai), path(assembly_fasta)

    output:
    tuple val(haplotype_id), path("${haplotype_id}.filtered.sorted.bam"), path("${haplotype_id}.filtered.sorted.bam.bai"), emit: bam
    tuple val(haplotype_id), path("${haplotype_id}_filtering_stats.txt"), emit: stats
    tuple val(haplotype_id), path("${haplotype_id}.pairs.gz"), emit: pairs

    script:
    """
    set -euo pipefail
    export LC_ALL=C

    MINQ=${params.hic_min_mapq ?: 30}
    TMPDIR="\${TMPDIR:-\$PWD}"

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
    samtools collate -@ ${task.cpus} -O -u ${bam} | \\
      pairtools parse \\
        --min-mapq \${MINQ} \\
        --walks-policy 5unique \\
        --max-inter-align-gap 30 \\
        --chroms-path chrom.sizes \\
        --output-stats ${haplotype_id}_parse_stats.txt \\
        --nproc-in  ${task.cpus} \\
        --nproc-out ${task.cpus} \\
        - \\
    | pairtools sort \\
        --nproc ${task.cpus} \\
        --tmpdir "\${TMPDIR}" \\
    | pairtools dedup \\
        --mark-dups \\
        --output-stats ${haplotype_id}_dedup_stats.txt \\
        --output-dups ${haplotype_id}.dups.pairs.gz \\
    | pairtools select '(pair_type == "UU")' \\
        --output ${haplotype_id}.pairsam.gz

    # -------------------------------------------------------------------------
    # 3) Create filtered .pairs.gz AND restore a filtered BAM, fast:
    #    pairtools split streams SAM to stdout; samtools consumes + sorts once.
    # -------------------------------------------------------------------------
    pairtools split \\
      --output-pairs ${haplotype_id}.pairs.gz \\
      --output-sam - \\
      --nproc-in  ${task.cpus} \\
      --nproc-out ${task.cpus} \\
      ${haplotype_id}.pairsam.gz \\
    | samtools view -@ ${task.cpus} -b - \\
    | samtools sort -@ ${task.cpus} -T "\${TMPDIR}/${haplotype_id}.tmp" -o ${haplotype_id}.filtered.sorted.bam -

    samtools index -@ ${task.cpus} ${haplotype_id}.filtered.sorted.bam

    # -------------------------------------------------------------------------
    # 4) Optional pair-level stats (pairtools stats) if available
    # -------------------------------------------------------------------------
    if pairtools stats --help >/dev/null 2>&1; then
      pairtools stats ${haplotype_id}.pairs.gz > ${haplotype_id}_pairs_stats.txt
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
    filtered_primary=\$(samtools view -c -F 2304 ${haplotype_id}.filtered.sorted.bam)

    # Avoid hard dependency on bc; awk is everywhere
    retained_pct=\$(awk -v f="\$filtered_primary" -v o="\$original_primary" 'BEGIN{ if(o>0) printf("%.2f", (f/o)*100); else print "NA"; }')

    {
      echo "================================================================================"
      echo "Hi-C BAM Filtering Statistics for ${haplotype_id}"
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
      cat ${haplotype_id}_parse_stats.txt
      echo
      echo "--------------------------------------------------------------------------------"
      echo "PAIRTOOLS DEDUP STATS"
      echo "--------------------------------------------------------------------------------"
      cat ${haplotype_id}_dedup_stats.txt
      echo
      if [[ "\${HAVE_PAIRTOOLS_STATS}" -eq 1 ]]; then
        echo "--------------------------------------------------------------------------------"
        echo "PAIRTOOLS PAIRS STATS (filtered pairs.gz)"
        echo "--------------------------------------------------------------------------------"
        cat ${haplotype_id}_pairs_stats.txt
        echo
      fi
      echo "--------------------------------------------------------------------------------"
      echo "FILTERED BAM (samtools flagstat)"
      echo "--------------------------------------------------------------------------------"
      samtools flagstat ${haplotype_id}.filtered.sorted.bam
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
    } > ${haplotype_id}_filtering_stats.txt

    # Cleanup
    rm -f ${haplotype_id}.pairsam.gz
    rm -f chrom.sizes
    """
    stub:
    """
    touch ${haplotype_id}.filtered.sorted.bam
    touch ${haplotype_id}.filtered.sorted.bam.bai
    touch ${haplotype_id}_filtering_stats.txt
    touch ${haplotype_id}.pairs.gz
    """
}
