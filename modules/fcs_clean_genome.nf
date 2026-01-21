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
  command -v fcs.py >/dev/null 2>&1 || { echo "ERROR: fcs.py not found in PATH"; exit 127; }

  # fcs.py clean genome reads from stdin in the canonical example
  cat ${assembly_fa} \\
    | fcs.py clean genome \\
        --action-report ${action_report} \\
        --output decontaminated.fasta \\
        --contam-fasta-out contaminants.fasta \\
    | tee clean_stdout.log
  """
}