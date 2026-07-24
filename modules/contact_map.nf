/*
========================================================================================
    CONTACT MAP MODULE (STANDALONE)
========================================================================================
    Repo location: modules/contact_map.nf

    Generates Hi-C contact maps from pre-filtered pairs.gz (no remapping).
    Input : tuple(meta, pairs_gz, assembly_fasta, stage)
    Output: cool / mcool / contact_map PNGs / stats  (per-haplotype, meta.id)
========================================================================================
*/

process CONTACT_MAP {
    tag "${meta.id}_${stage}"
    label 'hic_contact_map'

    publishDir "${params.outdir}/contact_maps", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(pairs_gz), path(assembly_fasta), val(stage)

    output:
    tuple val(meta), val(stage), path("${meta.id}_${stage}.cool"), emit: cool
    tuple val(meta), val(stage), path("${meta.id}_${stage}.mcool"), emit: mcool
    tuple val(meta), val(stage), path("${meta.id}_${stage}_*bp_contact_map.png"), emit: contact_maps
    tuple val(meta), val(stage), path("${meta.id}_${stage}_contact_stats.txt"), emit: stats

    script:
    def resolutions       = (params.hic_resolutions ?: "2500000,1000000,500000,250000,100000,50000,25000,10000").toString()
    def base_bin          = (params.hic_base_bin ?: "10000").toString()
    def plot_resolutions  = (params.hic_plot_resolutions ?: "1000000,500000,250000,100000").toString()
    def do_balance        = (params.run_hic_balance ?: false) as boolean
    def min_scaffold_size = params.scaffold_min_bp ?: 0  // 0 = include all scaffolds

    """
    set -euo pipefail

    # -------------------------------------------------------------------------
    # 1) Ensure reference index + chrom sizes
    # -------------------------------------------------------------------------
    if [ ! -s "${assembly_fasta}.fai" ]; then
        samtools faidx ${assembly_fasta}
    fi
    cut -f1,2 ${assembly_fasta}.fai > all_chrom.sizes

    # -------------------------------------------------------------------------
    # 2) Optional: Filter to major scaffolds only (for cleaner visualization)
    # -------------------------------------------------------------------------
    MIN_SIZE=${min_scaffold_size}

    if [[ \${MIN_SIZE} -gt 0 ]]; then
        echo "[CONTACT_MAP] Filtering to scaffolds >= \${MIN_SIZE} bp"

        # Create filtered chromosome sizes (major scaffolds only)
        awk -v min=\${MIN_SIZE} '\$2 >= min {print \$1"\\t"\$2}' all_chrom.sizes > chrom.sizes

        # Count filtered scaffolds
        TOTAL_SCAFFOLDS=\$(wc -l < all_chrom.sizes)
        MAJOR_SCAFFOLDS=\$(wc -l < chrom.sizes)
        echo "[CONTACT_MAP] Using \${MAJOR_SCAFFOLDS} of \${TOTAL_SCAFFOLDS} scaffolds (>= \${MIN_SIZE} bp)"

        # Build chromosome lookup for filtering pairs
        awk '{print \$1}' chrom.sizes > major_scaffolds.txt

        # Filter pairs.gz to only include major scaffolds
        zcat ${pairs_gz} \\
            | awk 'BEGIN{
                     while((getline < "major_scaffolds.txt") > 0) chr[\$1]=1
                   }
                   /^#/ {print; next}
                   chr[\$2] && chr[\$4] {print}' \\
            | bgzip -c > filtered.pairs.gz

        PAIRS_TO_USE="filtered.pairs.gz"
    else
        echo "[CONTACT_MAP] Using all scaffolds (no size filter)"
        cp all_chrom.sizes chrom.sizes
        PAIRS_TO_USE="${pairs_gz}"
    fi

    # -------------------------------------------------------------------------
    # 3) pairs.gz -> .cool at base resolution
    # -------------------------------------------------------------------------
    cooler cload pairs \\
        -c1 2 -p1 3 -c2 4 -p2 5 \\
        chrom.sizes:${base_bin} \\
        \${PAIRS_TO_USE} \\
        ${meta.id}_${stage}.cool

    # -------------------------------------------------------------------------
    # 4) .cool -> .mcool (multi-resolution), optionally balanced
    # -------------------------------------------------------------------------
    if ${do_balance}; then
        cooler zoomify \\
            --balance \\
            --resolutions ${resolutions} \\
            --out ${meta.id}_${stage}.mcool \\
            ${meta.id}_${stage}.cool
    else
        cooler zoomify \\
            --resolutions ${resolutions} \\
            --out ${meta.id}_${stage}.mcool \\
            ${meta.id}_${stage}.cool
    fi

    # -------------------------------------------------------------------------
    # 5) Generate contact map plots at specified resolutions
    # -------------------------------------------------------------------------
    IFS=',' read -r -a PLOT_RES <<< "${plot_resolutions}"
    for res in "\${PLOT_RES[@]}"; do
        hicPlotMatrix \\
            --matrix ${meta.id}_${stage}.mcool::/resolutions/\${res} \\
            --outFileName ${meta.id}_${stage}_\${res}bp_contact_map.png \\
            --title "${meta.id} (${stage}) Hi-C Contact Map (\${res} bp)" \\
            --log1p \\
            --dpi 200
    done

    # -------------------------------------------------------------------------
    # 6) Generate statistics report
    # -------------------------------------------------------------------------
    {
        echo "================================================================================"
        echo "Hi-C Contact Map Statistics for ${meta.id} (${stage})"
        echo "Generated: \$(date)"
        echo "================================================================================"
        echo
        echo "=== Input ==="
        echo "pairs.gz: ${pairs_gz}"
        echo "assembly_fasta: ${assembly_fasta}"
        echo
        if [[ \${MIN_SIZE} -gt 0 ]]; then
            echo "=== Scaffold Filtering ==="
            echo "Minimum scaffold size: \${MIN_SIZE} bp"
            echo "Total scaffolds in assembly: \${TOTAL_SCAFFOLDS}"
            echo "Major scaffolds used for visualization: \${MAJOR_SCAFFOLDS}"
            echo
            echo "Filtered scaffolds:"
            cat chrom.sizes
            echo
        fi
        echo "=== Pairs Summary ==="
        echo "Total pairs: \$(zcat \${PAIRS_TO_USE} | grep -v '^#' | wc -l)"
        echo
        echo "=== Cooler info (.cool) ==="
        cooler info ${meta.id}_${stage}.cool
        echo
        echo "=== Cooler tree (.mcool) ==="
        cooler tree ${meta.id}_${stage}.mcool
        echo
        echo "=== Resolutions ==="
        echo "Base bin size: ${base_bin}"
        echo "All resolutions: ${resolutions}"
        echo "Plot resolutions: ${plot_resolutions}"
        echo "Balancing: ${do_balance}"
    } > ${meta.id}_${stage}_contact_stats.txt
    """

    stub:
    """
    touch ${meta.id}_${stage}.cool
    touch ${meta.id}_${stage}.mcool
    touch ${meta.id}_${stage}_1000000bp_contact_map.png
    touch ${meta.id}_${stage}_contact_stats.txt
    """
}
