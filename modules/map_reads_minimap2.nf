process MAP_READS_MINIMAP2 {
  tag { meta?.id ?: "map_reads" }
  label 'mapping_qc' 

  input:
    path assembly_fa
    path reads
    val preset

  output:
    path "reads.bam",     emit: bam
    path "reads.bam.bai", emit: bai

  script:
  """
  set -euo pipefail

  minimap2 -t ${task.cpus} -ax ${preset} ${assembly_fa} ${reads} \
    | samtools sort -@ ${task.cpus} -o reads.bam

  samtools index -@ ${task.cpus} reads.bam
  """
}