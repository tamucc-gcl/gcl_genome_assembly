/*
========================================================================================
    HIFIASM MODULE
========================================================================================
    Genome assembly using HiFi (+ optional Hi-C) reads with hifiasm.
    Repo location: modules/hifiasm.nf

    Ploidy is meta-driven (Phase 2):
      - n_hap == 2 (diploid): phased hap1 + hap2 (Hi-C phasing if params.hifiasm_use_hic).
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
    def useHiC   = !haploid && params.hifiasm_use_hic && meta.hic
    def primary  = haploid || params.hifiasm_primary

    primary_flag  = primary ? '--primary' : ''
    dualscaf_flag = params.hifiasm_dual_scaffolding ? '--dual-scaf' : ''
    hic_opts      = useHiC ? "--h1 ${hic_r1} --h2 ${hic_r2}" : ''

    // Handle 'auto' parameters - omit flag entirely for auto behavior
    hgsize_opt   = params.hifiasm_haploid_genome_size   == 'auto' ? '' : "--hg-size ${params.hifiasm_haploid_genome_size}"
    homcov_opt   = params.hifiasm_homozygous_coverage   == 'auto' ? '' : "--hom-cov ${params.hifiasm_homozygous_coverage}"
    purgemax_opt = params.hifiasm_purge_max_coverage == 'auto' ? '' : "--purge-max ${params.hifiasm_purge_max_coverage}"

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
        -k ${params.hifiasm_kmer_length} \\
        -w ${params.hifiasm_minimizer_window} \\
        -f ${params.hifiasm_bloom_filter_bits} \\
        -D ${params.hifiasm_kmer_drop_factor} \\
        -N ${params.hifiasm_max_overlaps} \\
        -r ${params.hifiasm_correction_rounds} \\
        -z ${params.hifiasm_adapter_trim_bp} \\
        --max-kocc ${params.hifiasm_max_kmer_occurrence} \\
        ${hgsize_opt} \\
        -a ${params.hifiasm_cleaning_rounds} \\
        -m ${params.hifiasm_contig_bubble_bp} \\
        -p ${params.hifiasm_unitig_bubble_bp} \\
        -n ${params.hifiasm_tip_unitig_reads} \\
        -x ${params.hifiasm_max_overlap_drop_ratio} \\
        -y ${params.hifiasm_min_overlap_drop_ratio} \\
        -u ${params.hifiasm_post_join} \\
        ${homcov_opt} \\
        --lowQ ${params.hifiasm_low_quality_percent} \\
        --b-cov ${params.hifiasm_break_low_coverage} \\
        --h-cov ${params.hifiasm_break_high_coverage} \\
        --m-rate ${params.hifiasm_break_mismatch_rate} \\
        ${primary_flag} \\
        --ctg-n ${params.hifiasm_tip_contig_reads} \\
        -l ${params.hifiasm_purge_level} \\
        -s ${params.hifiasm_purge_similarity} \\
        -O ${params.hifiasm_purge_min_overlap} \\
        ${purgemax_opt} \\
        --n-hap ${params.hifiasm_nHaplotypes} \\
        ${dualscaf_flag} \\
        --scaf-gap ${params.hifiasm_scaffold_max_gap_bp} \\
        --telo-m ${telomere_motif} \\
        --telo-p ${params.hifiasm_telo_penalty} \\
        --telo-d ${params.hifiasm_telo_max_drop} \\
        --telo-s ${params.hifiasm_telo_min_score} \\
        --s-base ${params.hifiasm_hic_base_similarity} \\
        --n-weight ${params.hifiasm_hic_reweight_rounds} \\
        --n-perturb ${params.hifiasm_hic_perturb_rounds} \\
        --f-perturb ${params.hifiasm_hic_perturb_fraction} \\
        --l-msjoin ${params.hifiasm_misjoin_detect_bp} \\
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
