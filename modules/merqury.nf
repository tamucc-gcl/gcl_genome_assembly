/*
========================================================================================
    MERQURY MODULE (OPTIMIZED)
========================================================================================
    K-mer based assembly evaluation using pre-built meryl database
    - Accepts pre-computed meryl database (built once from HiFi reads)
    - Dramatically faster than rebuilding k-mer database each time
========================================================================================
*/

process MERQURY {
    tag "${sample_id}"
    label 'merqury'
    
    //publishDir "${params.outdir}/qc/assembly/merqury", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hap1_fasta), path(hap2_fasta), path(meryl_db)
    
    output:
    tuple val(sample_id), path("${sample_id}_merqury"), emit: results
    
    script:
    """
    mkdir -p ${sample_id}_merqury
    cd ${sample_id}_merqury
    
    # Run Merqury using pre-built k-mer database
    merqury.sh \\
        ../${meryl_db} \\
        ../${hap1_fasta} \\
        ../${hap2_fasta} \\
        ${sample_id}
    """
    
    stub:
    """
    mkdir -p ${sample_id}_merqury
    touch ${sample_id}_merqury/${sample_id}.qv
    touch ${sample_id}_merqury/${sample_id}.completeness.stats
    """
}