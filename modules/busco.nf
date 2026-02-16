/*
========================================================================================
    BUSCO MODULE
========================================================================================
    Assess genome completeness using BUSCO
    Note: Set params.busco_lineage in config (e.g., "actinopterygii_odb10")
========================================================================================
*/

process BUSCO {
    tag "${haplotype_id}"
    label 'busco'
    
    //publishDir "${params.outdir}/qc/assembly/busco", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly_fasta)
    path(busco_db)
    
    output:
    tuple val(haplotype_id), path("${haplotype_id}_busco"), emit: results
    
    script:
    def lineage = params.busco_lineage ?: "auto"
    """
    busco \\
        --in ${assembly_fasta} \\
        --out ${haplotype_id}_busco \\
        --mode genome \\
        --lineage_dataset ${lineage} \\
        --cpu ${task.cpus} \\
        --offline \\
        --download_path ${params.busco_downloads ?: '/path/to/busco_downloads'}
    """
    
    stub:
    """
    mkdir -p ${haplotype_id}_busco
    touch ${haplotype_id}_busco/short_summary.txt
    """
}