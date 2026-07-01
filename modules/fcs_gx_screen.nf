process FCS_GX_SCREEN {
  tag "${meta.id}"
  label 'fcs'

  publishDir "${params.outdir}/decontam/${stage}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> filename.startsWith('gx_out/') ? null : filename }

  input:
    tuple val(meta), path(assembly_fa), val(source_taxid), val(gxdb_dir), val(stage)

  output:
    tuple val(meta), path("*.fcs_gx_report.txt"), emit: action_report
    tuple val(meta), path("*.taxonomy.rpt"),      emit: taxonomy_report
    tuple val(meta), path("gx_out/fcs_gx_stdout.log"), emit: stdout_log

  stub:
  """
  mkdir -p gx_out
  touch ${meta.id}.fcs_gx_report.txt ${meta.id}.taxonomy.rpt gx_out/fcs_gx_stdout.log
  """

  script:
  """
  set -euo pipefail

  mkdir -p gx_out

  # Find the database prefix (e.g., test-only, all, etc.)
  DB_PREFIX=\$(ls ${gxdb_dir}/*.gxi | head -n1 | sed 's/\\.gxi\$//')

  if [ -z "\${DB_PREFIX}" ]; then
    echo "ERROR: No .gxi database file found in ${gxdb_dir}"
    ls -lh ${gxdb_dir}
    exit 1
  fi

  echo "Using GX database: \${DB_PREFIX}"
  echo "Processing haplotype: ${meta.id}"

  # Nextflow automatically wraps this in singularity exec
  /app/bin/run_gx \\
    --fasta ${assembly_fa} \\
    --out-dir gx_out \\
    --gx-db "\${DB_PREFIX}" \\
    --tax-id ${source_taxid} \\
    | tee gx_out/fcs_gx_stdout.log

  # Move output files to working directory for publishing
  # Add id to filenames for clarity
  mv gx_out/*.fcs_gx_report.txt ${meta.id}.fcs_gx_report.txt || true
  mv gx_out/*.taxonomy.rpt ${meta.id}.taxonomy.rpt || true
  """
}
