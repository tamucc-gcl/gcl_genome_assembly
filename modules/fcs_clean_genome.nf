process FCS_CLEAN_GENOME {
  tag "fcs_clean"
  label 'fcs' 

  input:
    path assembly_fa
    path action_report

  output:
    path "decontaminated.fasta", emit: decontaminated_fasta
    path "contaminants.fasta",   emit: contaminants_fasta
    path "clean_stdout.log",     emit: stdout_log

  script:
  """
  set -euo pipefail

  # Nextflow automatically wraps this in singularity exec
  /app/bin/action_report \\
    --fasta ${assembly_fa} \\
    --in ${action_report} \\
    --output decontaminated.fasta \\
    --contam-fasta-out contaminants.fasta \\
    | tee clean_stdout.log
  """
}