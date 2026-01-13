/*
========================================================================================
    HIFIASM MODULE
========================================================================================
    Phased genome assembly using HiFi and Hi-C reads with hifiasm
========================================================================================
*/

process HIFIASM {
    tag "${sample_id}"
    label 'hifiasm'
    
    publishDir "${params.outdir}/${sample_id}/assembly/hifiasm", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hifi_fastq), path(hic_r1), path(hic_r2)
    
    output:
    tuple val(sample_id), path("${sample_id}.hic.hap1.p_ctg.fasta.gz"), path("${sample_id}.hic.hap2.p_ctg.fasta.gz"), emit: assemblies
    tuple val(sample_id), path("${sample_id}.hifiasm.log"), emit: log
    tuple val(sample_id), path("${sample_id}.hic.hap1.p_ctg.gfa"), emit: gfa_hap1
    tuple val(sample_id), path("${sample_id}.hic.hap2.p_ctg.gfa"), emit: gfa_hap2
    
    script:
    """
    hifiasm \\
        -o ${sample_id} \\
        -t ${task.cpus} \\
        -l 3 \\
        -s 0.55 \\
        --n-hap 2 \\
        --h1 ${hic_r1} \\
        --h2 ${hic_r2} \\
        ${hifi_fastq} \\
        2>&1 | tee ${sample_id}.hifiasm.log
    
    # Convert GFA to FASTA
    gfatools gfa2fa ${sample_id}.hic.hap1.p_ctg.gfa > ${sample_id}.hic.hap1.p_ctg.fasta
    gfatools gfa2fa ${sample_id}.hic.hap2.p_ctg.gfa > ${sample_id}.hic.hap2.p_ctg.fasta
    
    # Compress FASTA files
    pigz -p ${task.cpus} -f ${sample_id}.hic.hap1.p_ctg.fasta
    pigz -p ${task.cpus} -f ${sample_id}.hic.hap2.p_ctg.fasta
    """
    
    stub:
    """
    touch ${sample_id}.hic.hap1.p_ctg.fasta.gz
    touch ${sample_id}.hic.hap2.p_ctg.fasta.gz
    touch ${sample_id}.hifiasm.log
    touch ${sample_id}.hic.hap1.p_ctg.gfa
    touch ${sample_id}.hic.hap2.p_ctg.gfa
    """
}