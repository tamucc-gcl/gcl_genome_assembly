/*
========================================================================================
    BUSCO MODULE
========================================================================================
    Repo location: modules/busco.nf

    Genome completeness per haplotype. Carries meta so ASSEMBLY_QC can regroup results
    per sample via groupKey(meta.sample, meta.n_hap).
    Note: set params.busco_lineage in config (e.g. "actinopterygii_odb10").
========================================================================================
*/

process BUSCO {
    tag "${meta.id}"
    label 'busco'

    //publishDir "${params.outdir}/qc/assembly/busco", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly_fasta)
    val(busco_lineage_map)

    output:
    tuple val(meta), path("${meta.id}_busco"), emit: results

    script:
    // Per-sample lineage, keyed by meta.taxid. Falls back to params.busco_lineage only if a
    // sample somehow lacks a resolved taxid (unreachable while taxid is required — see notes).
    def lineage = busco_lineage_map[meta.taxid?.toString()] ?: (params.busco_lineage ?: 'eukaryota_odb10')
    """
    busco \\
        --in ${assembly_fasta} \\
        --out ${meta.id}_busco \\
        --mode genome \\
        --lineage_dataset ${lineage} \\
        --cpu ${task.cpus} \\
        --offline \\
        --download_path ${params.busco_downloads ?: '/path/to/busco_downloads'}
    """

    stub:
    """
    mkdir -p ${meta.id}_busco
    touch ${meta.id}_busco/short_summary.txt
    """
}
