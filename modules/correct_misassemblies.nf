/*
========================================================================================
    MISASSEMBLY CORRECTION MODULE (INSPECTOR)
========================================================================================
    Identifies and corrects assembly errors using Inspector
    
    Workflow:
    1. inspector.py: Maps HiFi reads back to assembly and identifies structural errors
    2. inspector-correct.py: Breaks assembly at identified error positions
    
    Can be applied to:
    - Contig assemblies (post-HIFIASM)
    - Scaffolded assemblies (post-YaHS)
    
    Outputs:
    - Corrected assembly FASTA
    - Structural error BED file
    - Read-to-contig BAM alignment
    - Summary statistics
========================================================================================
*/

process CORRECT_MISASSEMBLIES {
    tag "${haplotype_id}"
    label 'misassembly_correction'
    
    publishDir "${params.outdir}/misassembly_correction/${stage}", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly_fasta), path(hifi_reads), val(stage), val(correction_params)
    
    output:
    tuple val(haplotype_id), path("${haplotype_id}_corrected.fa"), emit: corrected
    tuple val(haplotype_id), path("${haplotype_id}_read_to_contig.bam"), emit: bam
    tuple val(haplotype_id), path("${haplotype_id}_structural_error.bed"), emit: errors
    tuple val(haplotype_id), path("${haplotype_id}_summary_statistics.txt"), emit: stats
    tuple val(haplotype_id), path("${haplotype_id}_inspector"), emit: inspector_dir
    
    script:
    // Extract parameters for this stage
    def min_depth = correction_params.min_depth ?: ""
    def min_contig_length = correction_params.min_contig_length ?: 10000
    def min_contig_length_assemblyerror = correction_params.min_contig_length_assemblyerror ?: 1000000
    def min_assembly_error_size = correction_params.min_assembly_error_size ?: 50
    def max_assembly_error_size = correction_params.max_assembly_error_size ?: 4000000
    
    // Build optional arguments for inspector.py
    def depth_arg = min_depth ? "--min_depth ${min_depth}" : ""
    
    """
    set -euo pipefail
    
    # Create inspector output directory
    mkdir -p ${haplotype_id}_inspector
    cd ${haplotype_id}_inspector
    
    # Step 1: Identify assembly errors
    echo "[INSPECTOR] Identifying assembly errors for ${haplotype_id} (${stage})"
    echo "[INSPECTOR] Started: \$(date)"
    echo "[INSPECTOR] Parameters:"
    echo "  min_depth: ${min_depth ?: 'default (20% of average depth)'}"
    echo "  min_contig_length: ${min_contig_length}"
    echo "  min_contig_length_assemblyerror: ${min_contig_length_assemblyerror}"
    echo "  min_assembly_error_size: ${min_assembly_error_size}"
    echo "  max_assembly_error_size: ${max_assembly_error_size}"
    
    inspector.py \\
        --contig ../${assembly_fasta} \\
        --read ../${hifi_reads} \\
        --datatype hifi \\
        --thread ${task.cpus} \\
        --min_contig_length ${min_contig_length} \\
        --min_contig_length_assemblyerror ${min_contig_length_assemblyerror} \\
        --min_assembly_error_size ${min_assembly_error_size} \\
        --max_assembly_error_size ${max_assembly_error_size} \\
        ${depth_arg} \\
        --outpath ./
    
    echo "[INSPECTOR] Error identification complete: \$(date)"
    
    # Step 2: Break assembly at identified error positions
    echo "[INSPECTOR] Correcting assembly errors"
    echo "[INSPECTOR] Started: \$(date)"
    
    inspector-correct.py \\
        --inspector ./ \\
        --datatype pacbio-hifi \\
        --outpath ./ \\
        --thread ${task.cpus}
    
    echo "[INSPECTOR] Correction complete: \$(date)"
    
    # Move outputs to working directory with haplotype_id prefix
    cd ..
    
    # Corrected assembly
    if [ -f ${haplotype_id}_inspector/contig_corrected.fa ]; then
        mv ${haplotype_id}_inspector/contig_corrected.fa ${haplotype_id}_corrected.fa
    else
        echo "[WARNING] No corrections made - using original assembly"
        cp ${assembly_fasta} ${haplotype_id}_corrected.fa
    fi
    
    # Read-to-contig alignment BAM
    if [ -f ${haplotype_id}_inspector/read_to_contig.bam ]; then
        mv ${haplotype_id}_inspector/read_to_contig.bam ${haplotype_id}_read_to_contig.bam
    else
        touch ${haplotype_id}_read_to_contig.bam
    fi
    
    # Structural error BED
    if [ -f ${haplotype_id}_inspector/structural_error.bed ]; then
        mv ${haplotype_id}_inspector/structural_error.bed ${haplotype_id}_structural_error.bed
    else
        # Create empty BED if no errors found
        echo "# No structural errors detected" > ${haplotype_id}_structural_error.bed
    fi
    
    # Summary statistics
    if [ -f ${haplotype_id}_inspector/summary_statistics ]; then
        mv ${haplotype_id}_inspector/summary_statistics ${haplotype_id}_summary_statistics.txt
    else
        # Create basic summary if file doesn't exist
        cat > ${haplotype_id}_summary_statistics.txt <<EOF
Inspector Summary for ${haplotype_id} (${stage})
Generated: \$(date)

Assembly: ${assembly_fasta}
Reads: ${hifi_reads}

Status: Analysis complete
Check inspector directory for detailed results
EOF
    fi
    
    # Generate comprehensive report
    cat > ${haplotype_id}_inspector_report.txt <<EOF
================================================================================
Inspector Misassembly Correction Report
================================================================================
Haplotype: ${haplotype_id}
Stage: ${stage}
Generated: \$(date)

Input Assembly: ${assembly_fasta}
Input Reads: ${hifi_reads}

================================================================================
STRUCTURAL ERRORS DETECTED
================================================================================
\$(if [ -s ${haplotype_id}_structural_error.bed ]; then
    echo "Number of errors: \$(grep -v '^#' ${haplotype_id}_structural_error.bed | wc -l)"
    echo ""
    echo "Error locations:"
    grep -v '^#' ${haplotype_id}_structural_error.bed | head -20
    if [ \$(grep -v '^#' ${haplotype_id}_structural_error.bed | wc -l) -gt 20 ]; then
        echo "... (showing first 20 errors)"
    fi
else
    echo "No structural errors detected"
fi)

================================================================================
CORRECTION SUMMARY
================================================================================
\$(cat ${haplotype_id}_summary_statistics.txt)

================================================================================
OUTPUT FILES
================================================================================
Corrected assembly: ${haplotype_id}_corrected.fa
Structural errors: ${haplotype_id}_structural_error.bed
Read alignment: ${haplotype_id}_read_to_contig.bam
Summary stats: ${haplotype_id}_summary_statistics.txt
Full inspector output: ${haplotype_id}_inspector/

================================================================================
ASSEMBLY STATISTICS
================================================================================
Original assembly:
\$(if command -v seqkit &> /dev/null; then
    seqkit stats ${assembly_fasta}
else
    echo "Number of sequences: \$(grep -c '^>' ${assembly_fasta})"
    echo "Total length: \$(grep -v '^>' ${assembly_fasta} | tr -d '\\n' | wc -c)"
fi)

Corrected assembly:
\$(if command -v seqkit &> /dev/null; then
    seqkit stats ${haplotype_id}_corrected.fa
else
    echo "Number of sequences: \$(grep -c '^>' ${haplotype_id}_corrected.fa)"
    echo "Total length: \$(grep -v '^>' ${haplotype_id}_corrected.fa | tr -d '\\n' | wc -c)"
fi)

================================================================================
EOF
    
    # Add report to inspector directory
    cp ${haplotype_id}_inspector_report.txt ${haplotype_id}_inspector/
    """
    
    stub:
    """
    mkdir -p ${haplotype_id}_inspector
    touch ${haplotype_id}_corrected.fa
    touch ${haplotype_id}_read_to_contig.bam
    touch ${haplotype_id}_structural_error.bed
    touch ${haplotype_id}_summary_statistics.txt
    touch ${haplotype_id}_inspector_report.txt
    """
}