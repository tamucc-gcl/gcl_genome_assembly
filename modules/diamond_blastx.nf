process DIAMOND_BLASTX {
  tag "diamond_blastx"
  label 'diamond' 

  input:
    path assembly_fa
    val diamond_db      // Changed from 'path' to 'val' - don't stage the huge DB
    val max_target_seqs
    val evalue

  output:
    path "diamond_hits.tsv", emit: out_hits

  script:
  """
  set -euo pipefail

  diamond blastx \
    -d ${diamond_db} \
    -q ${assembly_fa} \
    -o diamond_hits.tsv \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
    --max-target-seqs ${(max_target_seqs ?: 1) as int} \
    --evalue ${evalue ?: 1e-25} \
    --threads ${task.cpus}
  """
}