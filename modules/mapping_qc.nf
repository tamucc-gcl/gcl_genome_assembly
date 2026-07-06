/*
========================================================================================
    MAPPING QC MODULE
========================================================================================
    Repo location: modules/mapping_qc.nf

    Maps the sample's QC reads back to each haplotype and computes coverage/mapping stats.
    Carries meta so ASSEMBLY_QC can regroup results per sample via groupKey.

    Phase 4a-iii: read-source-aware. `reads` is a single HiFi FASTQ for HiFi samples, or
    the Illumina R1+R2 pair for short-read samples; the minimap2 preset is chosen from meta
    (map-hifi for HiFi, sr for short-read paired-end, map-ont reserved for future long-read).
    `minimap2 -ax sr R1 R2` performs paired mapping.
========================================================================================
*/

process MAPPING_QC {
    tag "${meta.id}"
    label 'mapping_qc'

    //publishDir "${params.outdir}/qc/assembly/mapping", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly_fasta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_mapping_stats"), emit: results

    script:
    def preset = meta.hifi ? 'map-hifi' : ( meta.long_reads ? 'map-ont' : 'sr' )
    """
    mkdir -p ${meta.id}_mapping_stats

    # Index assembly
    minimap2 -d ${meta.id}.mmi ${assembly_fasta}

    # Map QC reads (preset by read type; `reads` = 1 HiFi FASTQ or R1 R2 for PE)
    minimap2 \\
        -ax ${preset} \\
        -t ${task.cpus} \\
        ${meta.id}.mmi \\
        ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam -

    # Index BAM
    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam

    # Calculate statistics
    samtools flagstat ${meta.id}.sorted.bam > ${meta.id}_mapping_stats/flagstat.txt
    samtools stats ${meta.id}.sorted.bam > ${meta.id}_mapping_stats/stats.txt
    samtools coverage ${meta.id}.sorted.bam > ${meta.id}_mapping_stats/coverage.txt

    # Calculate average depth
    awk 'NR==1{next} \$1!="*"{
        L = \$3 - \$2 + 1
        sum += \$7 * L
        len += L
        } END {
        print sum / len
        }' ${meta.id}_mapping_stats/coverage.txt > ${meta.id}_mapping_stats/avg_depth.txt
    """

    stub:
    """
    mkdir -p ${meta.id}_mapping_stats
    touch ${meta.id}_mapping_stats/flagstat.txt
    touch ${meta.id}_mapping_stats/stats.txt
    touch ${meta.id}_mapping_stats/coverage.txt
    touch ${meta.id}_mapping_stats/depth_summary.txt
    touch ${meta.id}_mapping_stats/avg_depth.txt
    """
}
