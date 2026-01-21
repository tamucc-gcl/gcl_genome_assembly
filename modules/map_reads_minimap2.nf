process MAP_READS_MINIMAP2 {
  tag { meta?.id ?: "map_reads" }
  label 'mapping_qc' 

  input:
    path assembly_fa
    path reads
    val preset
    val cpus

  output:
    path "reads.bam",     emit: bam
    path "reads.bam.bai", emit: bai

  // This minimap2 biocontainer *usually* includes samtools.
  // If your environment separates them, swap to a container image that has both.
  container "quay.io/biocontainers/minimap2:2.28--he4a0461_1"

  script:
  """
  set -euo pipefail

  minimap2 -t ${task.cpus} -ax ${preset} ${assembly_fa} ${reads} \
    | samtools sort -@ ${task.cpus} -o reads.bam

  samtools index -@ ${task.cpus} reads.bam
  """
}
