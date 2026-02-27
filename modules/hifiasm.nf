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
    
    publishDir "${params.outdir}/assembly/contig/hifiasm", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(hifi_fastq), path(hic_r1), path(hic_r2)
    
    output:
    tuple val(sample_id), path("${sample_id}.hap1.p_ctg.fasta"), path("${sample_id}.hap2.p_ctg.fasta"), emit: assemblies
    tuple val(sample_id), path("${sample_id}.hifiasm.log"), emit: log
    tuple val(sample_id), path("${sample_id}.hap1.p_ctg.gfa"), emit: gfa_hap1
    tuple val(sample_id), path("${sample_id}.hap2.p_ctg.gfa"), emit: gfa_hap2
    
    script:
    telomere_motif = params.telomere_motif ?: 'CCCTAA'  // Default to human telomeric repeat if not provided
    primary_flag = params.hifiasm_primary ? '--primary' : ''
    dualscaf_flag = params.hifiasm_dualScaf ? '--dual-scaf' : ''
    hic_opts = params.hifiasm_useHiC ? "--h1 ${hic_r1} --h2 ${hic_r2}" : ''

    // Determine source GFA naming pattern
    // HiFi-only + primary: *.p_ctg.gfa / *.a_ctg.gfa (no hap1/hap2)
    // HiFi-only default:   *.bp.hap1.p_ctg.gfa / *.bp.hap2.p_ctg.gfa
    // Hi-C (any):          *.hic.hap1.p_ctg.gfa / *.hic.hap2.p_ctg.gfa
    
    // Handle 'auto' parameters - omit flag entirely for auto behavior
    hgsize_opt = params.hifiasm_hgSize == 'auto' ? '' : "--hg-size ${params.hifiasm_hgSize}"
    homcov_opt = params.hifiasm_homCov == 'auto' ? '' : "--hom-cov ${params.hifiasm_homCov}"
    purgemax_opt = params.hifiasm_purgeMax == 'auto' ? '' : "--purge-max ${params.hifiasm_purgeMax}"


    use_primary_alt = params.hifiasm_primary && !params.hifiasm_useHiC
    prefix = params.hifiasm_useHiC ? 'hic.' : 'bp.'
    
    gfa1 = use_primary_alt ? "${sample_id}.p_ctg.gfa" : "${sample_id}.${prefix}hap1.p_ctg.gfa"
    gfa2 = use_primary_alt ? "${sample_id}.a_ctg.gfa" : "${sample_id}.${prefix}hap2.p_ctg.gfa"

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
        ${hgsize_opt} \\
        -a ${params.hifiasm_a} \\
        -m ${params.hifiasm_m} \\
        -p ${params.hifiasm_p} \\
        -n ${params.hifiasm_n} \\
        -x ${params.hifiasm_x} \\
        -y ${params.hifiasm_y} \\
        -u ${params.hifiasm_u} \\
        ${homcov_opt} \\
        --lowQ ${params.hifiasm_lowQ} \\
        --b-cov ${params.hifiasm_bCov} \\
        --h-cov ${params.hifiasm_hCov} \\
        --m-rate ${params.hifiasm_mRate} \\
        ${primary_flag} \\
        --ctg-n ${params.hifiasm_ctgN} \\
        -l ${params.hifiasm_l} \\
        -s ${params.hifiasm_s} \\
        -O ${params.hifiasm_O} \\
        ${purgemax_opt} \\
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
    gfatools gfa2fa ${gfa1} > ${sample_id}.hap1.p_ctg.fasta
    gfatools gfa2fa ${gfa2} > ${sample_id}.hap2.p_ctg.fasta

    # Copy GFAs to standardized names (preserves originals too)
    cp ${gfa1} ${sample_id}.hap1.p_ctg.gfa
    cp ${gfa2} ${sample_id}.hap2.p_ctg.gfa
    """
    
    stub:
    """
    touch ${sample_id}.hap1.p_ctg.fasta
    touch ${sample_id}.hap2.p_ctg.fasta
    touch ${sample_id}.hifiasm.log
    touch ${sample_id}.hap1.p_ctg.gfa
    touch ${sample_id}.hap2.p_ctg.gfa
    """
}