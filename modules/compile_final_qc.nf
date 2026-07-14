/*
========================================================================================
    COMPILE FINAL QC MODULE
========================================================================================
    Aggregates all QC metrics from across the pipeline into a single report
    
    Inputs:
    - Assembly QC summaries from all stages (contig, scaffold, gap-filled, etc.)
    - Hi-C BAM metrics from all checkpoints
    - Hi-C pairs metrics from all checkpoints
    
    Outputs:
    - Consolidated TSV report with all samples and stages
    - Optional HTML report for visualization
========================================================================================
*/

process COMPILE_FINAL_QC {
    tag "compile_final_qc"
    label 'summarize_assembly'
    
    publishDir "${params.outdir}/qc/assembly", mode: params.publish_dir_mode
    
    input:
    path(assembly_summaries)  // All assembly QC summary TSVs collected
    path(bam_metrics)         // All BAM metrics TSVs collected
    path(pairs_metrics)       // All pairs metrics TSVs collected
    path(compile_qc_script)   // R script for compiling QC data
    
    output:
    //path("final_qc_report.tsv"), emit: report
    //path("final_qc_report.html"), emit: html, optional: true
    //path("qc_inputs"), emit: inputs_dir
    path "assembly_qc_metrics.csv",  emit: metrics
    path "*.png",                     emit: plots
    
    script:
    """
    set -euo pipefail
    
    # Create organized input directory structure
    mkdir -p qc_inputs/assembly qc_inputs/bam qc_inputs/pairs
    
    # Copy assembly summaries (handle both single file and multiple files)
    for f in ${assembly_summaries}; do
        if [ -f "\$f" ]; then
            cp "\$f" qc_inputs/assembly/
        fi
    done
    
    # Copy BAM metrics
    for f in ${bam_metrics}; do
        if [ -f "\$f" ]; then
            cp "\$f" qc_inputs/bam/
        fi
    done
    
    # Copy pairs metrics
    for f in ${pairs_metrics}; do
        if [ -f "\$f" ]; then
            cp "\$f" qc_inputs/pairs/
        fi
    done
    
    # Create manifest of available inputs
    cat > qc_inputs/manifest.txt <<EOF
Pipeline Final QC Compilation
Generated: \$(date)

Assembly QC files:
\$(ls -1 qc_inputs/assembly/ 2>/dev/null || echo "  (none)")

BAM metrics files:
\$(ls -1 qc_inputs/bam/ 2>/dev/null || echo "  (none)")

Pairs metrics files:
\$(ls -1 qc_inputs/pairs/ 2>/dev/null || echo "  (none)")
EOF
    
    # Run R script to compile all QC data
    Rscript ${compile_qc_script} \\
        --assembly_dir qc_inputs/assembly \\
        --bam_dir qc_inputs/bam \\
        --pairs_dir qc_inputs/pairs \\
        --output_dir .
    """
    
    stub:
    """
    mkdir -p qc_inputs/assembly qc_inputs/bam qc_inputs/pairs
    touch assembly_qc_metrics.csv
    touch assembly_qc_metrics.png
    touch qc_inputs/manifest.txt
    """
}