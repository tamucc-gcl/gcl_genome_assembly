process FCS_GX_SCREEN {
  tag "fcs_gx_screen"
  label 'fcs' 

  input:
    path assembly_fa
    val source_taxid
    path gxdb_dir

  output:
    path "gx_out/*.fcs_gx_report.txt", emit: action_report
    path "gx_out/*.taxonomy.rpt",      emit: taxonomy_report
    path "gx_out/fcs_gx_stdout.log",   emit: stdout_log

  cpus { cpus ?: 32 }

  script:
  """
  command -v fcs.py >/dev/null 2>&1 || { echo "ERROR: fcs.py not found in PATH"; exit 127; }

  mkdir -p gx_out

  export GX_NUM_CORES=${task.cpus}

  # If your input is gzipped fasta, fcs.py can handle it (as shown in docs),
  # but plain fasta is fine too.
  fcs.py screen genome \
    --fasta ${assembly_fa} \
    --out-dir gx_out \
    --gx-db ${gxdb_dir} \
    --tax-id ${source_taxid} \
    | tee gx_out/fcs_gx_stdout.log
  """
}
