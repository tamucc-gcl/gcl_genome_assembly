/*
========================================================================================
    FINALIZE ASSEMBLY MODULE
========================================================================================
    Takes the final assembly from the pipeline and publishes it to a clean
    output directory with a standardized filename.

    This module serves as the single collection point for "the final genome"
    regardless of which optional pipeline steps were run upstream:
      - Gap filling only (default)
      - Gap filling + teloclip extension (if params.run_teloclip_extend)

    The input channel is ch_final_assembly, which is set in main.nf to either
    TELOCLIP_EXTEND.out.extended_assembly or GAP_FILLING.out.filled_assembly.

    Input:
    - tuple(haplotype_id, assembly_fasta)

    Output:
    - tuple(haplotype_id, final_fasta)
========================================================================================
*/

process FINALIZE_ASSEMBLY {
    tag "${haplotype_id}"
    label 'finalize_assembly'

    publishDir "${params.outdir}/assembly/final", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(assembly_fasta, stageAs: 'input/*')

    output:
    tuple val(haplotype_id), path("${haplotype_id}.fasta"), emit: assembly

    script:
    """
    cp ${assembly_fasta.name} ${haplotype_id}.fasta
    """

    stub:
    """
    touch ${haplotype_id}.fasta
    """
}