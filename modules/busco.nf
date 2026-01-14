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
    
    publishDir "${params.outdir}/qc/assembly/busco", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly_fasta)
    
    output:
    tuple val(haplotype_id), path("${haplotype_id}_busco"), emit: results
    
    script:
    def lineage = params.busco_lineage ?: "auto"
    """
        # Decompress assemblies if needed
    if [[ ${assembly_fasta} =~ \\.gz\$ ]]; then
        pigz -dc ${assembly_fasta} > hap.fasta
        hap_input="hap.fasta"
    else
        hap_input="${hap1_fasta}"
    fi

    busco \\
        --in ${hap_input} \\
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