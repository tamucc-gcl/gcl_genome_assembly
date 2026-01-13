/*
========================================================================================
    MAPPING QC MODULE
========================================================================================
    Map HiFi reads back to assembly and calculate coverage/mapping statistics
========================================================================================
*/

process MAPPING_QC {
    tag "${haplotype_id}"
    label 'mapping_qc'
    
    publishDir "${params.outdir}/${haplotype_id}/qc/assembly/mapping", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id), path(assembly_fasta), path(hifi_fastq)
    
    output:
    tuple val(haplotype_id), path("${haplotype_id}_mapping_stats"), emit: results
    
    script:
    """
    mkdir -p ${haplotype_id}_mapping_stats
    
    # Index assembly
    minimap2 -d ${haplotype_id}.mmi ${assembly_fasta}
    
    # Map HiFi reads
    minimap2 \\
        -ax map-hifi \\
        -t ${task.cpus} \\
        ${haplotype_id}.mmi \\
        ${hifi_fastq} \\
        | samtools sort -@ ${task.cpus} -o ${haplotype_id}.sorted.bam -
    
    # Index BAM
    samtools index -@ ${task.cpus} ${haplotype_id}.sorted.bam
    
    # Calculate statistics
    samtools flagstat ${haplotype_id}.sorted.bam > ${haplotype_id}_mapping_stats/flagstat.txt
    samtools stats ${haplotype_id}.sorted.bam > ${haplotype_id}_mapping_stats/stats.txt
    samtools coverage ${haplotype_id}.sorted.bam > ${haplotype_id}_mapping_stats/coverage.txt
    
    # Calculate average depth
    awk 'NR==1{next} \$1!="*"{
        L = \$3 - \$2 + 1
        sum += \$7 * L
        len += L
        } END {
        print sum / len
        }' ${haplotype_id}_mapping_stats/coverage.txt > ${haplotype_id}_mapping_stats/avg_depth.txt
    """
    
    stub:
    """
    mkdir -p ${haplotype_id}_mapping_stats
    touch ${haplotype_id}_mapping_stats/flagstat.txt
    touch ${haplotype_id}_mapping_stats/stats.txt
    touch ${haplotype_id}_mapping_stats/coverage.txt
    touch ${haplotype_id}_mapping_stats/depth_summary.txt
    touch ${haplotype_id}_mapping_stats/avg_depth.txt
    """
}