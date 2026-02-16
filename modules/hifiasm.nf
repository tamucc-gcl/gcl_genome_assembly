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
    
    publishDir "${params.outdir}/contig/hifiasm", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hifi_fastq), path(hic_r1), path(hic_r2)
    
    output:
    tuple val(sample_id), path("${sample_id}.hic.hap1.p_ctg.fasta"), path("${sample_id}.hic.hap2.p_ctg.fasta"), emit: assemblies
    tuple val(sample_id), path("${sample_id}.hifiasm.log"), emit: log
    tuple val(sample_id), path("${sample_id}.hic.hap1.p_ctg.gfa"), emit: gfa_hap1
    tuple val(sample_id), path("${sample_id}.hic.hap2.p_ctg.gfa"), emit: gfa_hap2
    
    script:
    telomere_motif = params.telomere_motif ?: 'CCCTAA'  // Default to human telomeric repeat if not provided
    primary_flag = params.hifiasm_primary ? '--primary' : ''
    dualscaf_flag = params.hifiasm_dualScaf ? '--dual-scaf' : ''
    hic_opts = params.hifiasm_useHiC ? "--h1 ${hic_r1} --h2 ${hic_r2}" : ''
    """
    hifiasm \\
        -o ${sample_id} \\
        -t ${task.cpus} \\
        -k ${params.hifiasm_k} \\
        -w ${params.hifiasm_w} \\
        -f ${params.hifiasm_f} \\
        -D ${params.hifiasm_D} \\
        -N ${params.hifiasm_N} \\
        -r ${params.hifiasm_r} \\
        -z ${params.hifiasm_z} \\
        --max-kocc ${params.hifiasm_maxKOCC} \\
        --hg-size ${params.hifiasm_hgSize} \\
        -a ${params.hifiasm_a} \\
        -m ${params.hifiasm_m} \\
        -p ${params.hifiasm_p} \\
        -n ${params.hifiasm_n} \\
        -x ${params.hifiasm_x} \\
        -y ${params.hifiasm_y} \\
        -u ${params.hifiasm_u} \\
        --hom-cov ${params.hifiasm_homCov} \\
        --lowQ ${params.hifiasm_lowQ} \\
        --b-cov ${params.hifiasm_bCov} \\
        --h-cov ${params.hifiasm_hCov} \\
        --m-rate ${params.hifiasm_mRate} \\
        ${primary_flag} \\
        --ctg-n ${params.hifiasm_ctgN} \\
        -l ${params.hifiasm_l} \\
        -s ${params.hifiasm_s} \\
        -O ${params.hifiasm_O} \\
        --purge-max ${params.hifiasm_purgeMax} \\
        --n-hap ${params.hifiasm_nHaplotypes} \\
        ${dualscaf_flag} \\
        --scaf-gap ${params.hifiasm_scafGap} \\
        --telo-m ${telomere_motif} \\
        --telo-p ${params.hifiasm_teloP} \\
        --telo-d ${params.hifiasm_teloD} \\
        --telo-s ${params.hifiasm_teloS} \\
        --s-base ${params.hifiasm_sBase} \\
        --n-weight ${params.hifiasm_nWeight} \\
        --n-perturb ${params.hifiasm_nPerturb} \\
        --f-perturb ${params.hifiasm_fPerturb} \\
        --l-msjoin ${params.hifiasm_lMSjoin} \\
        ${hic_opts} \\
        ${hifi_fastq} \\
        2>&1 | tee ${sample_id}.hifiasm.log
    
    # Convert GFA to FASTA
    gfatools gfa2fa ${sample_id}.hic.hap1.p_ctg.gfa > ${sample_id}.hic.hap1.p_ctg.fasta
    gfatools gfa2fa ${sample_id}.hic.hap2.p_ctg.gfa > ${sample_id}.hic.hap2.p_ctg.fasta
    """
    
    stub:
    """
    touch ${sample_id}.hic.hap1.p_ctg.fasta
    touch ${sample_id}.hic.hap2.p_ctg.fasta
    touch ${sample_id}.hifiasm.log
    touch ${sample_id}.hic.hap1.p_ctg.gfa
    touch ${sample_id}.hic.hap2.p_ctg.gfa
    """
}