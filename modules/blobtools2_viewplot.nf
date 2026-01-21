process BLOBTOOLS_VIEWPLOT {
  tag "blobtools_viewplot"
  label 'blobtools' 

  input:
    path blobdir

  output:
    path "blobtools_outputs", emit: out_dir

  script:
  """
  set -euo pipefail

  mkdir -p blobtools_outputs

  # Human-readable evidence
  blobtools view ${blobdir} -o blobtools_outputs
  blobtools plot ${blobdir} -o blobtools_outputs

  # Also try to export a machine-readable table if supported by your blobtools2 build.
  # If not supported, this step will be skipped (but view/plot outputs remain).
  if blobtools --help 2>/dev/null | grep -qE 'table\\b'; then
    blobtools table ${blobdir} > blobtools_outputs/blobtable.tsv || true
  fi
  """
}
