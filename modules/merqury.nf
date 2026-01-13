/*
========================================================================================
    MERQURY MODULE
========================================================================================
    K-mer based assembly evaluation using HiFi reads
========================================================================================
*/

process MERQURY {
    tag "${sample_id}"
    label 'merqury'
    
    publishDir "${params.outdir}/qc/assembly/merqury", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hap1_fasta), path(hap2_fasta), path(hifi_fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}_merqury"), emit: results
    
    script:
    """
    # Decompress assemblies if needed
    if [[ ${hap1_fasta} =~ \\.gz\$ ]]; then
        pigz -dc ${hap1_fasta} > hap1.fasta
        hap1_input="hap1.fasta"
    else
        hap1_input="${hap1_fasta}"
    fi
    
    if [[ ${hap2_fasta} =~ \\.gz\$ ]]; then
        pigz -dc ${hap2_fasta} > hap2.fasta
        hap2_input="hap2.fasta"
    else
        hap2_input="${hap2_fasta}"
    fi

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