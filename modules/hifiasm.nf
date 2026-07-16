/*
========================================================================================
    HIFIASM MODULE
========================================================================================
    Genome assembly using HiFi (+ optional Hi-C) reads with hifiasm.
    Repo location: modules/hifiasm.nf

    Ploidy is meta-driven (Phase 2):
      - n_hap == 2 (diploid): phased hap1 + hap2 (Hi-C phasing if params.hifiasm_useHiC).
        Emits two FASTAs: <sample>.hap1.p_ctg.fasta, <sample>.hap2.p_ctg.fasta.
      - n_hap == 1 (haploid/collapsed): forces --primary and DISABLES Hi-C phasing
        (Hi-C, if present, is still used downstream for scaffolding the single assembly).
        Emits one FASTA: <sample>.primary.p_ctg.fasta (from the p_ctg graph).

    Output `assemblies` is therefore variable-arity: a single FASTA (haploid) or two
    (diploid). main.nf zips forkHaplotypeMeta(meta) against it (coercing to a list), so
    haploid yields one (primary, fasta) and diploid two (hapN, fasta) tuples.
========================================================================================
*/

process HIFIASM {
    tag "${meta.sample}"
    label 'hifiasm'

    publishDir "${params.outdir}/assembly/contig/hifiasm", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hifi_fastq), path(hic_r1), path(hic_r2)

    output:
    tuple val(meta), path("${meta.sample}.{hap1,hap2,primary}.p_ctg.fasta"), emit: assemblies
    tuple val(meta), path("${meta.sample}.hifiasm.log"),                     emit: log
    tuple val(meta), path("${meta.sample}.{hap1,hap2,primary}.p_ctg.gfa"),   emit: gfa

    script:
    telomere_motif = params.telomere_motif ?: 'CCCTAA'  // Default to human telomeric repeat if not provided

    // ---- ploidy-driven mode (Phase 2) ----
    // Haploid -> collapsed primary assembly: --primary, no Hi-C phasing.
    // Hi-C phasing also requires Hi-C to actually be present (HiFi-only rows -> off).
    def haploid  = (meta.n_hap == 1)
    def useHiC   = !haploid && params.hifiasm_useHiC && meta.hic
    def primary  = haploid || params.hifiasm_primary

    primary_flag  = primary ? '--primary' : ''
    dualscaf_flag = params.hifiasm_dualScaf ? '--dual-scaf' : ''
    hic_opts      = useHiC ? "--h1 ${hic_r1} --h2 ${hic_r2}" : ''

    // Handle 'auto' parameters - omit flag entirely for auto behavior
    hgsize_opt   = params.hifiasm_hgSize   == 'auto' ? '' : "--hg-size ${params.hifiasm_hgSize}"
    homcov_opt   = params.hifiasm_homCov   == 'auto' ? '' : "--hom-cov ${params.hifiasm_homCov}"
    purgemax_opt = params.hifiasm_purgeMax == 'auto' ? '' : "--purge-max ${params.hifiasm_purgeMax}"

    // Source GFA naming:
    //   primary (no Hi-C):  *.p_ctg.gfa / *.a_ctg.gfa
    //   HiFi-only default:  *.bp.hap1.p_ctg.gfa / *.bp.hap2.p_ctg.gfa
    //   Hi-C phased:        *.hic.hap1.p_ctg.gfa / *.hic.hap2.p_ctg.gfa
    use_primary_alt = primary && !useHiC
    prefix = useHiC ? 'hic.' : 'bp.'
    gfa1 = use_primary_alt ? "${meta.sample}.p_ctg.gfa" : "${meta.sample}.${prefix}hap1.p_ctg.gfa"
    gfa2 = use_primary_alt ? "${meta.sample}.a_ctg.gfa" : "${meta.sample}.${prefix}hap2.p_ctg.gfa"

    // Haploid: single primary output from gfa1 (== <sample>.p_ctg.gfa).
    // Diploid: two haplotype outputs from gfa1/gfa2.
    def conversion_cmds
    if (haploid) {
        conversion_cmds = "gfatools gfa2fa ${gfa1} > ${meta.sample}.primary.p_ctg.fasta\n    cp ${gfa1} ${meta.sample}.primary.p_ctg.gfa"
    } else {
        conversion_cmds = "gfatools gfa2fa ${gfa1} > ${meta.sample}.hap1.p_ctg.fasta\n    gfatools gfa2fa ${gfa2} > ${meta.sample}.hap2.p_ctg.fasta\n    cp ${gfa1} ${meta.sample}.hap1.p_ctg.gfa\n    cp ${gfa2} ${meta.sample}.hap2.p_ctg.gfa"
    }

    """
    hifiasm \\
        -o ${meta.sample} \\
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
        2>&1 | tee ${meta.sample}.hifiasm.log

    # Convert GFA to FASTA (+ copy GFAs to standardized names)
    ${conversion_cmds}
    """

    stub:
    def haploid = (meta.n_hap == 1)
    if (haploid)
        """
        touch ${meta.sample}.primary.p_ctg.fasta
        touch ${meta.sample}.primary.p_ctg.gfa
        touch ${meta.sample}.hifiasm.log
        """
    else
        """
        touch ${meta.sample}.hap1.p_ctg.fasta
        touch ${meta.sample}.hap2.p_ctg.fasta
        touch ${meta.sample}.hap1.p_ctg.gfa
        touch ${meta.sample}.hap2.p_ctg.gfa
        touch ${meta.sample}.hifiasm.log
        """
}
