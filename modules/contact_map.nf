/*
========================================================================================
    CONTACT MAP MODULE (STANDALONE)
========================================================================================
    Generates Hi-C contact maps from pre-filtered pairs.gz files
    
    Input:
    - haplotype_id: Sample/haplotype identifier
    - pairs_gz: Filtered pairs.gz from FILTER_HIC_BAM
    - assembly_fasta: Reference assembly (contigs, scaffolds, or gap-filled)
    - stage: Label for output organization (e.g., "contig", "scaffold", "final")
    
    Output:
    - cool: Single-resolution cooler file
    - mcool: Multi-resolution cooler file
    - contact_maps: PNG visualizations at specified resolutions
    - stats: Contact map statistics
    
    Key features:
    - No remapping required - uses pre-computed pairs
    - Optional filtering to major scaffolds for cleaner visualization
    - Configurable resolutions and balancing
========================================================================================
*/

process CONTACT_MAP {
    tag "${haplotype_id}_${stage}"
    label 'hic_contact_map'

    publishDir "${params.outdir}/contact_maps", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(pairs_gz), path(assembly_fasta), val(stage)

    output:
    tuple val(haplotype_id), val(stage), path("${haplotype_id}_${stage}.cool"), emit: cool
    tuple val(haplotype_id), val(stage), path("${haplotype_id}_${stage}.mcool"), emit: mcool
    tuple val(haplotype_id), val(stage), path("${haplotype_id}_${stage}_*bp_contact_map.png"), emit: contact_maps
    tuple val(haplotype_id), val(stage), path("${haplotype_id}_${stage}_contact_stats.txt"), emit: stats

    script:
    def resolutions       = (params.hic_resolutions ?: "1000000,500000,100000,50000,10000").toString()
    def base_bin          = (params.hic_base_bin ?: "10000").toString()
    def plot_resolutions  = (params.hic_plot_resolutions ?: "1000000,500000,100000").toString()
    def do_balance        = (params.hic_balance ?: false) as boolean
    def min_scaffold_size = params.scaffold_min_size ?: 0  // 0 = include all scaffolds
    
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
        ${haplotype_id}_${stage}.cool

    # -------------------------------------------------------------------------
    # 4) .cool -> .mcool (multi-resolution), optionally balanced
    # -------------------------------------------------------------------------
    if ${do_balance}; then
        cooler zoomify \\
            --balance \\
            --resolutions ${resolutions} \\
            --out ${haplotype_id}_${stage}.mcool \\
            ${haplotype_id}_${stage}.cool
    else
        cooler zoomify \\
            --resolutions ${resolutions} \\
            --out ${haplotype_id}_${stage}.mcool \\
            ${haplotype_id}_${stage}.cool
    fi

    # -------------------------------------------------------------------------
    # 5) Generate contact map plots at specified resolutions
    # -------------------------------------------------------------------------
    IFS=',' read -r -a PLOT_RES <<< "${plot_resolutions}"
    for res in "\${PLOT_RES[@]}"; do
        hicPlotMatrix \\
            --matrix ${haplotype_id}_${stage}.mcool::/resolutions/\${res} \\
            --outFileName ${haplotype_id}_${stage}_\${res}bp_contact_map.png \\
            --title "${haplotype_id} (${stage}) Hi-C Contact Map (\${res} bp)" \\
            --log1p \\
            --dpi 200
    done

    # -------------------------------------------------------------------------
    # 6) Generate statistics report
    # -------------------------------------------------------------------------
    {
        echo "================================================================================"
        echo "Hi-C Contact Map Statistics for ${haplotype_id} (${stage})"
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
        cooler info ${haplotype_id}_${stage}.cool
        echo
        echo "=== Cooler tree (.mcool) ==="
        cooler tree ${haplotype_id}_${stage}.mcool
        echo
        echo "=== Resolutions ==="
        echo "Base bin size: ${base_bin}"
        echo "All resolutions: ${resolutions}"
        echo "Plot resolutions: ${plot_resolutions}"
        echo "Balancing: ${do_balance}"
    } > ${haplotype_id}_${stage}_contact_stats.txt
    """

    stub:
    """
    touch ${haplotype_id}_${stage}.cool
    touch ${haplotype_id}_${stage}.mcool
    touch ${haplotype_id}_${stage}_1000000bp_contact_map.png
    touch ${haplotype_id}_${stage}_contact_stats.txt
    """
}