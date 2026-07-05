/*
========================================================================================
    MERQURY MODULE (OPTIMIZED)
========================================================================================
    K-mer based assembly evaluation using a pre-built meryl database.
    Variable haplotype count (Phase 2): merqury.sh accepts one (haploid) or two
    (diploid) assemblies. Caller (assembly_qc.nf) passes the ordered FASTA list.
========================================================================================
*/

process MERQURY {
    tag "${sample_id}"
    label 'merqury'

    //publishDir "${params.outdir}/qc/assembly/merqury", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(fastas), path(meryl_db)

    output:
    tuple val(sample_id), path("${sample_id}_merqury"), emit: results

    script:
    // Assemblies are staged in the task dir; merqury runs from a subdir, so prefix with ../
    def asm_args = (fastas instanceof List ? fastas : [fastas]).collect { "../${it}" }.join(' ')
    """
    mkdir -p ${sample_id}_merqury
    cd ${sample_id}_merqury

    # Run Merqury using pre-built k-mer database (1 assembly = haploid, 2 = diploid)
    merqury.sh \\
        ../${meryl_db} \\
        ${asm_args} \\
        ${sample_id}
    """

    stub:
    """
    mkdir -p ${sample_id}_merqury
    touch ${sample_id}_merqury/${sample_id}.qv
    touch ${sample_id}_merqury/${sample_id}.completeness.stats
    """
}
