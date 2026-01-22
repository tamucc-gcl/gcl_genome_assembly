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
  command -v run_gx.py >/dev/null 2>&1 || { echo "ERROR: run_gx.py not found in PATH"; exit 127; }

  # run_gx.py clean genome reads from stdin in the canonical example
  cat ${assembly_fa} \\
    | run_gx.py clean genome \\
        --action-report ${action_report} \\
        --output decontaminated.fasta \\
        --contam-fasta-out contaminants.fasta \\
    | tee clean_stdout.log
  """
}