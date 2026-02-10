/*
========================================================================================
    COMBINE ASSEMBLY QC MODULE
========================================================================================
    Aggregates all assembly QC results and generates summary plots using R
========================================================================================
*/

process COMBINE_ASSEMBLY_QC {
    tag "${sample_id}:${qc_label}"
    label 'summarize_assembly'
    
    publishDir "${params.outdir}/qc/assembly/${qc_label}/", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), 
          val(qc_label),
          path(quast_results), 
          path(merqury_results),
          val(haplotype_ids_busco),
          path(busco_results),
          val(haplotype_ids_mapping),
          path(mapping_results)
    
    output:
    tuple val(sample_id), val(qc_label), path("${sample_id}_${qc_label}_qc_summary.tsv"), emit: summary
    
    script:
    """
    # Create directory structure for all QC inputs
    mkdir -p ${sample_id}_qc_inputs/quast
    mkdir -p ${sample_id}_qc_inputs/merqury
    mkdir -p ${sample_id}_qc_inputs/busco
    mkdir -p ${sample_id}_qc_inputs/mapping
    
    # Copy QUAST results
    cp -r ${quast_results}/* ${sample_id}_qc_inputs/quast/
    
    # Copy only the MERQURY stats files (not the huge .meryl databases!)
    cp ${merqury_results}/*.qv ${sample_id}_qc_inputs/merqury/ 2>/dev/null || true
    cp ${merqury_results}/*.completeness.stats ${sample_id}_qc_inputs/merqury/ 2>/dev/null || true
    cp ${merqury_results}/*.png ${sample_id}_qc_inputs/merqury/ 2>/dev/null || true
    cp ${merqury_results}/*.hist ${sample_id}_qc_inputs/merqury/ 2>/dev/null || true
    
    # Copy BUSCO results - only the summary files
    for busco_dir in ${busco_results}; do
        busco_name=\$(basename \$busco_dir)
        mkdir -p ${sample_id}_qc_inputs/busco/\$busco_name
        cp \$busco_dir/short_summary.*.json ${sample_id}_qc_inputs/busco/\$busco_name/ 2>/dev/null || true
        cp \$busco_dir/short_summary.*.txt ${sample_id}_qc_inputs/busco/\$busco_name/ 2>/dev/null || true
    done
    
    # Copy mapping results - only the stats files
    for mapping_dir in ${mapping_results}; do
        mapping_name=\$(basename \$mapping_dir)
        mkdir -p ${sample_id}_qc_inputs/mapping/\$mapping_name
        cp \$mapping_dir/*.txt ${sample_id}_qc_inputs/mapping/\$mapping_name/
    done
    
    # Create a manifest file listing all inputs
    cat > ${sample_id}_qc_inputs/manifest.txt <<EOF
SAMPLE_ID: ${sample_id}
QC_LABEL: ${qc_label}
QUAST_DIR: quast
MERQURY_DIR: merqury
BUSCO_HAPLOTYPES: ${haplotype_ids_busco.join(',')}
BUSCO_DIRS: \$(ls -d ${sample_id}_qc_inputs/busco/* | xargs -n1 basename | paste -sd,)
MAPPING_HAPLOTYPES: ${haplotype_ids_mapping.join(',')}
MAPPING_DIRS: \$(ls -d ${sample_id}_qc_inputs/mapping/* | xargs -n1 basename | paste -sd,)
EOF
    
    
    # Call R script to join everything within a sample together
    Rscript ${projectDir}/r_scripts/combine_individual_assembly_qc.R \\
        --input_dir ${sample_id}_qc_inputs \\
        --output_dir . \\
        --sample_id ${sample_id} \\
        --qc_label ${qc_label}
    """
    
    stub:
    """
    mkdir -p ${sample_id}_qc_inputs
    touch ${sample_id}_${qc_label}_qc_summary.tsv
    """
}