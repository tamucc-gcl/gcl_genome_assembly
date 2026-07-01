/*
========================================================================================
    MISASSEMBLY CORRECTION MODULE (INSPECTOR)
========================================================================================
    Identifies and corrects assembly errors using Inspector
    Repo location: modules/correct_misassemblies.nf

    Workflow:
    1. inspector.py: Maps HiFi reads back to assembly and identifies structural errors
    2. inspector-correct.py: Breaks assembly at identified error positions

    Reused for contig (post-HIFIASM) and scaffold (post-YaHS) stages via the `stage` val.

    Outputs:
    - Corrected assembly FASTA, structural-error BED, read-to-contig BAM, summary stats
========================================================================================
*/

process CORRECT_MISASSEMBLIES {
    tag "${meta.id}"
    label 'misassembly_correction'

    publishDir "${params.outdir}/assembly/${stage}/misassembly_correction", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly_fasta), path(hifi_reads), val(stage), val(correction_params)

    output:
    tuple val(meta), path("${meta.id}_corrected.fasta"), emit: corrected
    tuple val(meta), path("${meta.id}_read_to_contig.bam"), emit: bam
    tuple val(meta), path("${meta.id}_structural_error.bed"), emit: structural_errors
    tuple val(meta), path("${meta.id}_small_scale_error.bed"), emit: small_scale_errors
    tuple val(meta), path("${meta.id}_summary_statistics.txt"), emit: stats
    tuple val(meta), path("${meta.id}_inspector_report.txt"), emit: inspector_report
    //tuple val(meta), path("${meta.id}_inspector"), emit: inspector_dir

    script:
    // Extract parameters for this stage
    def skip_base = correction_params.skip_baseerror ? "--skip_baseerror" : ""
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
    mkdir -p ${meta.id}_inspector
    cd ${meta.id}_inspector

    # Step 1: Identify assembly errors
    echo "[INSPECTOR] Identifying assembly errors for ${meta.id} (${stage})"
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
        ${skip_base} \\
        --datatype pacbio-hifi \\
        --outpath ./ \\
        --thread ${task.cpus}

    echo "[INSPECTOR] Correction complete: \$(date)"

    # Move outputs to working directory with id prefix
    cd ..

    # Corrected assembly
    if [ -f ${meta.id}_inspector/contig_corrected.fa ]; then
        mv ${meta.id}_inspector/contig_corrected.fa ${meta.id}_corrected.fasta
    else
        echo "[WARNING] No corrections made - using original assembly"
        cp ${assembly_fasta} ${meta.id}_corrected.fasta
    fi

    # Read-to-contig alignment BAM
    if [ -f ${meta.id}_inspector/read_to_contig.bam ]; then
        mv ${meta.id}_inspector/read_to_contig.bam ${meta.id}_read_to_contig.bam
    else
        touch ${meta.id}_read_to_contig.bam
    fi

    # Structural error BED
    if [ -f ${meta.id}_inspector/structural_error.bed ]; then
        mv ${meta.id}_inspector/structural_error.bed ${meta.id}_structural_error.bed
    else
        # Create empty BED if no errors found
        echo "# No structural errors detected" > ${meta.id}_structural_error.bed
    fi

    # Small scale error BED
    if [ -f ${meta.id}_inspector/small_scale_error.bed ]; then
        mv ${meta.id}_inspector/small_scale_error.bed ${meta.id}_small_scale_error.bed
    else
        # Create empty BED if no errors found
        echo "# No small scale errors detected" > ${meta.id}_small_scale_error.bed
    fi

    # Summary statistics
    if [ -f ${meta.id}_inspector/summary_statistics ]; then
        mv ${meta.id}_inspector/summary_statistics ${meta.id}_summary_statistics.txt
    else
        # Create basic summary if file doesn't exist
        cat > ${meta.id}_summary_statistics.txt <<EOF
Inspector Summary for ${meta.id} (${stage})
Generated: \$(date)

Assembly: ${assembly_fasta}
Reads: ${hifi_reads}

Status: Analysis complete
Check inspector directory for detailed results
EOF
    fi

    # Generate comprehensive report
    cat > ${meta.id}_inspector_report.txt <<EOF
================================================================================
Inspector Misassembly Correction Report
================================================================================
Haplotype: ${meta.id}
Stage: ${stage}
Generated: \$(date)

Input Assembly: ${assembly_fasta}
Input Reads: ${hifi_reads}

================================================================================
STRUCTURAL ERRORS DETECTED
================================================================================
\$(if [ -s ${meta.id}_structural_error.bed ]; then
    echo "Number of errors: \$(grep -v '^#' ${meta.id}_structural_error.bed | wc -l)"
    echo ""
    echo "Error locations:"
    grep -v '^#' ${meta.id}_structural_error.bed | head -20
    if [ \$(grep -v '^#' ${meta.id}_structural_error.bed | wc -l) -gt 20 ]; then
        echo "... (showing first 20 errors)"
    fi
else
    echo "No structural errors detected"
fi)

================================================================================
CORRECTION SUMMARY
================================================================================
\$(cat ${meta.id}_summary_statistics.txt)

================================================================================
OUTPUT FILES
================================================================================
Corrected assembly: ${meta.id}_corrected.fasta
Structural errors: ${meta.id}_structural_error.bed
Read alignment: ${meta.id}_read_to_contig.bam
Summary stats: ${meta.id}_summary_statistics.txt

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
    seqkit stats ${meta.id}_corrected.fasta
else
    echo "Number of sequences: \$(grep -c '^>' ${meta.id}_corrected.fasta)"
    echo "Total length: \$(grep -v '^>' ${meta.id}_corrected.fasta | tr -d '\\n' | wc -c)"
fi)

================================================================================
EOF

    """

    stub:
    """
    mkdir -p ${meta.id}_inspector
    touch ${meta.id}_corrected.fasta
    touch ${meta.id}_read_to_contig.bam
    touch ${meta.id}_structural_error.bed
    touch ${meta.id}_small_scale_error.bed
    touch ${meta.id}_summary_statistics.txt
    touch ${meta.id}_inspector_report.txt
    """
}
