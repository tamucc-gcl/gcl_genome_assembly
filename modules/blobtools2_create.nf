process BLOBTOOLS_CREATE {
  tag "blobtools_create"
  label 'blobtools' 

  input:
    path assembly_fa
    path diamond_hits
    path bam
    path taxdump_dir
    val min_len

  output:
    path "blobdir", emit: out_blobdir

  script:
  """
  set -euo pipefail

  mkdir -p blobdir

  blobtools create \
    --fasta ${assembly_fa} \
    --bam ${bam} \
    --hits ${diamond_hits} \
    --taxdump ${taxdump_dir} \
    --min-length ${(min_len ?: 1000) as int} \
    blobdir
  """
}
