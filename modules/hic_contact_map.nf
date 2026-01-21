/*
========================================================================================
    HI-C CONTACT MAP MODULE (FIXED + FASTER + MORE EXPLICIT)
========================================================================================
    BAM -> pairs.gz -> cool -> mcool (+ optional balance) -> plots + stats

    Key improvements:
    - Fixes broken pipeline for qc_label != "raw"
    - Explicitly selects UU pairs (better “valid pair” semantics)
    - Uses TMPDIR for pairtools sort temp
    - Optional balancing (zoomify --balance can be very expensive)
    - Configurable plotting resolutions
    - Avoids re-running faidx when .fai exists
========================================================================================
*/

process HIC_CONTACT_MAP {
    tag "${haplotype_id}_${qc_label}"
    label 'hic_contact_map'

    publishDir "${params.outdir}/qc/hic_mapping/${qc_label}/${haplotype_id}/contact_maps",
        mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(bam), path(bai), path(assembly_fasta), val(qc_label)

    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.cool"), emit: cool
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.pairs.gz"), emit: pairs
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.mcool"), emit: mcool
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_*bp_contact_map.png"), emit: contact_maps
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_contact_stats.txt"), emit: stats

    script:
    // Resolutions to build in mcool
    def resolutions         = (params.hic_resolutions ?: "1000000,500000,100000,50000,10000").toString()
    // Base bin size used for initial .cool (must be one of the resolutions; typically the smallest)
    def base_bin            = (params.hic_base_bin ?: "10000").toString()
    // Resolutions to plot
    def plot_resolutions    = (params.hic_plot_resolutions ?: "1000000,500000,100000").toString()
    // Make balancing optional (default: false for speed; set true when you need it)
    def do_balance          = (params.hic_balance ?: false) as boolean
    // MapQ: raw maps often need strict; filtered can be lower because you already filtered upstream
    def min_mapq_raw        = (params.hic_min_mapq_raw ?: 30) as int
    def min_mapq_filtered   = (params.hic_min_mapq_filtered ?: 1) as int

    """
    set -euo pipefail
    export LC_ALL=C
    TMPDIR="\${TMPDIR:-\$PWD}"

    # -------------------------------------------------------------------------
    # 1) chrom sizes (avoid reindexing if already present)
    # -------------------------------------------------------------------------
    if [[ ! -s "${assembly_fasta}.fai" ]]; then
      samtools faidx ${assembly_fasta}
    fi
    cut -f1,2 ${assembly_fasta}.fai > chrom.sizes

    # -------------------------------------------------------------------------
    # 2) BAM -> pairs.gz
    #    Always end with a command that actually writes pairs output.
    #    We keep only UU pairs to make maps reflect high-confidence contacts.
    # -------------------------------------------------------------------------
    if [[ "${qc_label}" == "raw" ]]; then
      MINQ=${min_mapq_raw}

      samtools collate -@ ${task.cpus} -O -u ${bam} | \\
        pairtools parse \\
            --min-mapq \${MINQ} \\
            --walks-policy 5unique \\
            --max-inter-align-gap 30 \\
            --chroms-path chrom.sizes \\
            --output-stats ${haplotype_id}_${qc_label}_parse_stats.txt \\
            --nproc-in  ${task.cpus} \\
            --nproc-out ${task.cpus} \\
            - \\
      | pairtools select '(pair_type == "UU")' \\
      | pairtools sort --nproc ${task.cpus} --tmpdir "\${TMPDIR}" \\
      | pairtools dedup --mark-dups \\
          --output-stats ${haplotype_id}_${qc_label}_dedup_stats.txt \\
      | pairtools split --output-pairs ${haplotype_id}_${qc_label}.pairs.gz
    else
      MINQ=${min_mapq_filtered}

      samtools collate -@ ${task.cpus} -O -u ${bam} | \\
        pairtools parse \\
          --min-mapq \${MINQ} \\
          --walks-policy all \\
          --chroms-path chrom.sizes \\
          --output-stats ${haplotype_id}_${qc_label}_parse_stats.txt \\
          --nproc-in  ${task.cpus} \\
          --nproc-out ${task.cpus} \\
          - \\
      | pairtools select '(pair_type == "UU")' \\
      | pairtools sort --nproc ${task.cpus} --tmpdir "\${TMPDIR}" \\
      | pairtools split --output-pairs ${haplotype_id}_${qc_label}.pairs.gz
    fi

    # -------------------------------------------------------------------------
    # 3) pairs.gz -> .cool at base resolution
    # -------------------------------------------------------------------------
    cooler cload pairs \\
      -c1 2 -p1 3 -c2 4 -p2 5 \\
      chrom.sizes:${base_bin} \\
      ${haplotype_id}_${qc_label}.pairs.gz \\
      ${haplotype_id}_${qc_label}.cool

    # -------------------------------------------------------------------------
    # 4) .cool -> .mcool (multi-resolution), optionally balanced
    # -------------------------------------------------------------------------
    if ${do_balance}; then
      cooler zoomify \\
        --balance \\
        --resolutions ${resolutions} \\
        --out ${haplotype_id}_${qc_label}.mcool \\
        ${haplotype_id}_${qc_label}.cool
    else
      cooler zoomify \\
        --resolutions ${resolutions} \\
        --out ${haplotype_id}_${qc_label}.mcool \\
        ${haplotype_id}_${qc_label}.cool
    fi

    # -------------------------------------------------------------------------
    # 5) Plots (only for resolutions requested)
    # -------------------------------------------------------------------------
    IFS=',' read -r -a PLOT_RES <<< "${plot_resolutions}"
    for res in "\${PLOT_RES[@]}"; do
      hicPlotMatrix \\
        --matrix ${haplotype_id}_${qc_label}.mcool::/resolutions/\${res} \\
        --outFileName ${haplotype_id}_${qc_label}_\${res}bp_contact_map.png \\
        --title "${haplotype_id} (${qc_label}) Hi-C Contact Map (\${res} bp)" \\
        --log1p \\
        --dpi 200
    done

    # -------------------------------------------------------------------------
    # 6) Stats
    # -------------------------------------------------------------------------
    {
      echo "================================================================================"
      echo "Hi-C Contact Map Statistics for ${haplotype_id} (${qc_label})"
      echo "Generated: \$(date)"
      echo "================================================================================"
      echo
      echo "=== Input ==="
      echo "bam: ${bam}"
      echo "assembly_fasta: ${assembly_fasta}"
      echo
      echo "=== Pairtools parse stats ==="
      cat ${haplotype_id}_${qc_label}_parse_stats.txt
      echo
      if [[ "${qc_label}" == "raw" ]]; then
        echo "=== Pairtools dedup stats ==="
        cat ${haplotype_id}_${qc_label}_dedup_stats.txt
        echo
      fi
      echo "=== Cooler info (.cool) ==="
      cooler info ${haplotype_id}_${qc_label}.cool
      echo
      echo "=== Cooler tree (.mcool) ==="
      cooler tree ${haplotype_id}_${qc_label}.mcool
    } > ${haplotype_id}_${qc_label}_contact_stats.txt
    """
}


/*
========================================================================================
    HI-C CONTACT MAP MODULE (PAIRS INPUT)
========================================================================================
    pairs.gz -> cool -> mcool (+ optional balance) -> plots + stats

    This variant allows reusing an existing pairs.gz (e.g., lifted to scaffold coordinates)
    without remapping or re-parsing a BAM.

    Note: The plotting/stats steps are intentionally kept consistent with HIC_CONTACT_MAP.
========================================================================================
*/

process HIC_CONTACT_MAP_FROM_PAIRS {
    tag "${haplotype_id}_${qc_label}"
    label 'hic_contact_map'

    publishDir "${params.outdir}/qc/hic_mapping/${qc_label}/${haplotype_id}/contact_maps", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(pairs_gz), path(assembly_fasta), val(qc_label)

    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.cool"), emit: cool
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.mcool"), emit: mcool
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_*bp_contact_map.png"), emit: contact_maps
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_contact_stats.txt"), emit: stats

    script:
    """
    set -euo pipefail

    # -------------------------------------------------------------------------
    # 0) Inputs / parameters
    # -------------------------------------------------------------------------
    TMPDIR="\${TMPDIR:-\$PWD}"
    base_bin=${params.hic_base_bin ?: 1000}
    zoom_bins="${params.hic_zoom_bins ?: '1000,5000,10000,25000,50000,100000'}"
    do_balance=${params.hic_balance ?: false}
    plot_bins="${params.hic_plot_bins ?: '10000,50000,100000'}"

    # -------------------------------------------------------------------------
    # 1) Ensure reference index + chrom sizes
    # -------------------------------------------------------------------------
    if [ ! -s "${assembly_fasta}.fai" ]; then
      samtools faidx ${assembly_fasta}
    fi
    cut -f1,2 ${assembly_fasta}.fai > chrom.sizes

    # -------------------------------------------------------------------------
    # 2) pairs.gz -> .cool at base resolution
    # -------------------------------------------------------------------------
    cooler cload pairs \
      -c1 2 -p1 3 -c2 4 -p2 5 \
      chrom.sizes:${base_bin} \
      ${pairs_gz} \
      ${haplotype_id}_${qc_label}.cool

    # -------------------------------------------------------------------------
    # 3) .cool -> .mcool (multi-resolution), optionally balanced
    # -------------------------------------------------------------------------
    if [[ "\${do_balance}" == "true" ]]; then
      cooler zoomify --balance \
        --resolutions ${zoom_bins} \
        --nproc ${task.cpus} \
        ${haplotype_id}_${qc_label}.cool \
        -o ${haplotype_id}_${qc_label}.mcool
    else
      cooler zoomify \
        --resolutions ${zoom_bins} \
        --nproc ${task.cpus} \
        ${haplotype_id}_${qc_label}.cool \
        -o ${haplotype_id}_${qc_label}.mcool
    fi

    # -------------------------------------------------------------------------
    # 4) Plot contact maps at selected resolutions
    # -------------------------------------------------------------------------
    for BIN in \$(echo ${plot_bins} | tr ',' ' '); do
      cooler show \
        --dpi 200 \
        --out ${haplotype_id}_${qc_label}_\${BIN}bp_contact_map.png \
        ${haplotype_id}_${qc_label}.mcool::/resolutions/\${BIN}
    done

    # -------------------------------------------------------------------------
    # 5) Basic contact map stats
    # -------------------------------------------------------------------------
    {
      echo "haplotype_id\tqc_label\tcool_file\tmcool_file"
      echo "${haplotype_id}\t${qc_label}\t${haplotype_id}_${qc_label}.cool\t${haplotype_id}_${qc_label}.mcool"
    } > ${haplotype_id}_${qc_label}_contact_stats.txt
    """

    stub:
    """
    touch ${haplotype_id}_${qc_label}.cool
    touch ${haplotype_id}_${qc_label}.mcool
    touch ${haplotype_id}_${qc_label}_10000bp_contact_map.png
    touch ${haplotype_id}_${qc_label}_contact_stats.txt
    """
}
