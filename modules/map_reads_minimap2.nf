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
    path "versions.tsv", emit: versions

  script:
  """
  set -euo pipefail

  minimap2 -t ${task.cpus} -ax ${preset} ${assembly_fa} ${reads} \
    | samtools sort -@ ${task.cpus} -o reads.bam

  samtools index -@ ${task.cpus} reads.bam

  printf 'minimap2\t%s\n' "$(minimap2 --version 2>&1)" > versions.tsv
  printf 'samtools\t%s\n' "$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')" >> versions.tsv
  """
}