process FCS_CLEAN_GENOME {
  tag "${meta.id}"
  label 'fcs'

  publishDir "${params.outdir}/assembly/${stage}/decontam", mode: params.publish_dir_mode

  input:
    tuple val(meta), path(assembly_fa), path(action_report), val(stage)

  output:
    tuple val(meta), path("${meta.id}.decontaminated.fasta"), emit: decontaminated_fasta
    tuple val(meta), path("${meta.id}.contaminants.fasta"),   emit: contaminants_fasta

  stub:
  """
  touch ${meta.id}.decontaminated.fasta ${meta.id}.contaminants.fasta
  """

  script:
  """
  set -euo pipefail

  /app/bin/gx clean-genome \\
    --action-report ${action_report} \\
    --input ${assembly_fa} \\
    --output ${meta.id}.decontaminated.fasta \\
    --contam-fasta-out ${meta.id}.contaminants.fasta
  """
}
