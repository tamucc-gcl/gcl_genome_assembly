/*
========================================================================================
    SUMMARY REPORT MODULE
========================================================================================
    Generates a comprehensive HTML summary report with:
    1. Visual table: rows=samples, cols=hap1 snail, hap2 snail, hap1 vs hap2 dotplot
    2. Summary QC metrics table (assembly stats, Hi-C mapping, pairs stats)
    3. Assembly QC plots (BUSCO, N50/contiguity, QV, k-mer completeness)
    4. Links to all relevant output files
========================================================================================
*/

process SUMMARY_REPORT {
    tag "summary_report"
    label 'summary_report'
    
    conda "${projectDir}/environments/report.yaml"
    
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode
    
    input:
    // Snail plots: tuple(haplotype_id, qc_label, snail_svg)
    path(snail_plots)
    // Dotplots: tuple(hap1_id, hap2_id, dotplot_png)
    path(dotplots)
    // Contact map images: tuple(haplotype_id, stage, contact_png)
    path(contact_map_images)
    // Assembly QC summaries: tuple(sample_id, qc_label, summary_tsv)
    path(assembly_summaries)
    // BAM metrics: tuple(haplotype_id, checkpoint, metrics_tsv)
    path(bam_metrics)
    // Pairs metrics: tuple(haplotype_id, checkpoint, metrics_tsv)
    path(pairs_metrics)
    // Final assemblies
    path(final_assemblies)
    // QUAST final report directory
    path(quast_report)
    
    output:
    path("pipeline_summary_report.html"), emit: html_report
    path("pipeline_summary_report.md"), emit: md_report
    path("summary_data/"), emit: data_dir
    
    script:
    """
    set -euo pipefail
    
    mkdir -p summary_data/snail_plots
    mkdir -p summary_data/dotplots
    mkdir -p summary_data/contact_maps
    mkdir -p summary_data/qc_metrics
    mkdir -p summary_data/assemblies
    
    # Copy snail plots
    for f in ${snail_plots}; do
        [ -f "\$f" ] && cp "\$f" summary_data/snail_plots/
    done
    
    # Copy dotplots
    for f in ${dotplots}; do
        [ -f "\$f" ] && cp "\$f" summary_data/dotplots/
    done
    
    # Copy contact map images
    for f in ${contact_map_images}; do
        [ -f "\$f" ] && cp "\$f" summary_data/contact_maps/
    done
    
    # Copy QC metrics
    for f in ${assembly_summaries}; do
        [ -f "\$f" ] && cp "\$f" summary_data/qc_metrics/
    done
    for f in ${bam_metrics}; do
        [ -f "\$f" ] && cp "\$f" summary_data/qc_metrics/
    done
    for f in ${pairs_metrics}; do
        [ -f "\$f" ] && cp "\$f" summary_data/qc_metrics/
    done
    
    # Generate the report
    Rscript ${projectDir}/r_scripts/generate_summary_report.R \
        --snail_dir summary_data/snail_plots \
        --dotplot_dir summary_data/dotplots \
        --contact_map_dir summary_data/contact_maps \
        --assembly_qc_dir summary_data/qc_metrics \
        --output_dir . \
        --outdir_base "${params.outdir}"
    """
    
    stub:
    """
    mkdir -p summary_data
    touch pipeline_summary_report.html
    touch pipeline_summary_report.md
    """
}