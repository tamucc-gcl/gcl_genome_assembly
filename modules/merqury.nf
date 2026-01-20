/*
========================================================================================
    MERQURY MODULE
========================================================================================
    K-mer based assembly evaluation using HiFi reads
========================================================================================
*/

process MERQURY {
    tag "${sample_id}_${qc_label}"
    label 'merqury'
    
    //publishDir "${params.outdir}/qc/assembly/${qc_label}/merqury", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hap1_fasta), path(hap2_fasta), path(hifi_fastq), val(qc_label)
    
    output:
    tuple val(sample_id), path("${sample_id}_merqury"), val(qc_label), emit: results
    
    script:
    """
    mkdir -p ${sample_id}_merqury
    cd ${sample_id}_merqury
    
    # Build k-mer database from HiFi reads
    meryl count \\
        k=21 \\
        threads=${task.cpus} \\
        memory=${task.memory.toGiga()} \\
        ../${hifi_fastq} \\
        output reads.meryl
    
    # Run Merqury
    merqury.sh \\
        reads.meryl \\
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