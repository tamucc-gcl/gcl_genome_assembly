process FCS_ADAPTOR {
  tag "${meta.id}"

  publishDir "${params.outdir}/${stage}/adaptor", mode: params.publish_dir_mode

  input:
    tuple val(meta), path(assembly_fa), val(mode), val(engine), val(stage)

  output:
    tuple val(meta), path("${meta.id}.cleaned.fasta"), emit: cleaned_fasta
    tuple val(meta), path("fcsadaptor"),               emit: out_dir

  stub:
  """
  mkdir -p fcsadaptor
  touch ${meta.id}.cleaned.fasta
  """

  script:
  """
  command -v run_fcsadaptor.sh >/dev/null 2>&1 || { echo "ERROR: run_fcsadaptor.sh not found in PATH"; exit 127; }

  mkdir -p fcsadaptor/input fcsadaptor/output
  cp ${assembly_fa} fcsadaptor/input/assembly.fa

  # For fish genomes: treat as eukaryote (no --prok flag)
  if [ "${mode}" = "prok" ]; then
    run_fcsadaptor.sh --fasta-input fcsadaptor/input/assembly.fa --output-dir fcsadaptor/output --prok --container-engine ${engine}
  else
    run_fcsadaptor.sh --fasta-input fcsadaptor/input/assembly.fa --output-dir fcsadaptor/output --container-engine ${engine}
  fi

  # normalize expected output name
  # (FCS-adaptor produces a cleaned fasta in its output dir)
  ls -1 fcsadaptor/output/*.fa* | head -n 1 > /tmp/_cleaned_path.txt
  CLEANED=\$(cat /tmp/_cleaned_path.txt)
  cp "\$CLEANED" ${meta.id}.cleaned.fasta
  """
}
