/*
========================================================================================
    COMBINE HI-C MAPPING QC MODULE
========================================================================================
    Aggregates all Hi-C mapping QC results and generates summary report
    Can process either raw or filtered BAM QC results
========================================================================================
*/

process COMBINE_HIC_MAPPING_QC {
    tag "${haplotype_id}_${qc_label}"
    label 'combine_hic_qc'
    
    publishDir "${params.outdir}/qc/hic_mapping/${qc_label}/${haplotype_id}/summary", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id),
          val(qc_label),
          path(contact_stats),
          path(pair_stats),
          path(coverage_stats),
          path(contact_map_plots),
          path(pair_plots),
          path(coverage_plots)
    
    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_hic_mapping_report.html"), emit: report
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_hic_mapping_summary.txt"), emit: summary
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}_hic_qc_inputs"), emit: inputs
    
    script:
    """
    # Create directory structure for all QC inputs
    mkdir -p ${haplotype_id}_${qc_label}_hic_qc_inputs/contact_maps
    mkdir -p ${haplotype_id}_${qc_label}_hic_qc_inputs/pair_stats
    mkdir -p ${haplotype_id}_${qc_label}_hic_qc_inputs/coverage
    
    # Copy all inputs to organized directory structure
    cp ${contact_stats} ${haplotype_id}_${qc_label}_hic_qc_inputs/contact_maps/
    cp ${pair_stats} ${haplotype_id}_${qc_label}_hic_qc_inputs/pair_stats/
    cp ${coverage_stats} ${haplotype_id}_${qc_label}_hic_qc_inputs/coverage/
    
    # Copy plots
    cp ${contact_map_plots} ${haplotype_id}_${qc_label}_hic_qc_inputs/contact_maps/ 2>/dev/null || true
    cp ${pair_plots} ${haplotype_id}_${qc_label}_hic_qc_inputs/pair_stats/ 2>/dev/null || true
    cp ${coverage_plots} ${haplotype_id}_${qc_label}_hic_qc_inputs/coverage/ 2>/dev/null || true
    
    # Create comprehensive summary
    cat > ${haplotype_id}_${qc_label}_hic_mapping_summary.txt <<EOF
================================================================================
Hi-C Mapping QC Summary for ${haplotype_id} (${qc_label})
================================================================================
Generated: \$(date)

BAM Type: ${qc_label}

================================================================================
CONTACT MAP STATISTICS
================================================================================
\$(cat ${contact_stats})

================================================================================
PAIR STATISTICS
================================================================================
\$(cat ${pair_stats})

================================================================================
COVERAGE STATISTICS
================================================================================
\$(cat ${coverage_stats})

================================================================================
QUALITY ASSESSMENT
================================================================================
EOF

    # Add quality assessment based on statistics
    # Extract key metrics and evaluate quality
    trans_cis=\$(grep "Trans/Cis ratio" ${pair_stats} | awk '{print \$NF}')
    valid_pairs=\$(grep "mapped (" ${pair_stats} | head -1 | awk '{print \$(NF-1)}')
    
    cat >> ${haplotype_id}_${qc_label}_hic_mapping_summary.txt <<EOF

Key Quality Metrics:
- BAM Type: ${qc_label}
- Valid pairs percentage: \${valid_pairs}
- Trans/Cis ratio: \${trans_cis}

Quality Indicators:
- Trans/Cis ratio < 0.2: Excellent phasing
- Trans/Cis ratio 0.2-0.5: Good phasing
- Trans/Cis ratio > 0.5: Poor phasing, possible misassemblies

Notes:
- Low trans/cis ratio indicates good haplotype separation
- Uniform coverage suggests even Hi-C library preparation
- High mapping rates indicate good assembly contiguity
- Raw BAM includes all reads; Filtered BAM includes only valid pairs

================================================================================
EOF

    # Placeholder for HTML report (to be implemented with R Markdown or similar)
    cat > ${haplotype_id}_${qc_label}_hic_mapping_report.html <<EOF
<!DOCTYPE html>
<html>
<head>
    <title>Hi-C Mapping QC Report - ${haplotype_id} (${qc_label})</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; }
        h2 { color: #3498db; border-bottom: 2px solid #3498db; }
        pre { background-color: #f4f4f4; padding: 10px; border-radius: 5px; }
        .metric { background-color: #e8f4f8; padding: 10px; margin: 10px 0; border-radius: 5px; }
        .badge { 
            display: inline-block; 
            padding: 5px 10px; 
            border-radius: 3px; 
            font-weight: bold; 
            margin-left: 10px;
        }
        .raw { background-color: #ffebcc; color: #cc6600; }
        .filtered { background-color: #ccffcc; color: #006600; }
    </style>
</head>
<body>
    <h1>Hi-C Mapping QC Report</h1>
    <h2>Sample: ${haplotype_id} <span class="badge ${qc_label}">${qc_label}</span></h2>
    <p>Generated: \$(date)</p>
    
    <h2>Summary</h2>
    <div class="metric">
        <p><strong>BAM Type:</strong> ${qc_label}</p>
        <p><strong>Valid Pairs:</strong> \${valid_pairs}</p>
        <p><strong>Trans/Cis Ratio:</strong> \${trans_cis}</p>
    </div>
    
    <h2>Contact Maps</h2>
    <p>See PNG files in contact_maps directory</p>
    
    <h2>Detailed Statistics</h2>
    <pre>\$(cat ${haplotype_id}_${qc_label}_hic_mapping_summary.txt)</pre>
    
</body>
</html>
EOF

    # Future: Call R script to generate comprehensive HTML report
    # Rscript \${projectDir}/bin/combine_hic_mapping_qc.R \\
    #     --input_dir ${haplotype_id}_${qc_label}_hic_qc_inputs \\
    #     --output_dir . \\
    #     --haplotype_id ${haplotype_id} \\
    #     --qc_label ${qc_label}
    """
    
    stub:
    """
    mkdir -p ${haplotype_id}_${qc_label}_hic_qc_inputs
    touch ${haplotype_id}_${qc_label}_hic_mapping_report.html
    touch ${haplotype_id}_${qc_label}_hic_mapping_summary.txt
    """
}