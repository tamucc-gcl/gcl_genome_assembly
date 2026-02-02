process FCS_ADAPTOR {
  tag "${haplotype_id}"  // ← ADDED: Tag with haplotype_id

  publishDir "${params.outdir}/${stage}/adaptor", mode: params.publish_dir_mode

  input:
    tuple val(haplotype_id), path(assembly_fa), val(mode), val(engine), val(stage)  // ← CHANGED: Tuple input

  output:
    tuple val(haplotype_id), path("${haplotype_id}.cleaned.fasta"), emit: cleaned_fasta  // ← CHANGED: Tuple output
    tuple val(haplotype_id), path("fcsadaptor"),                    emit: out_dir

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
  cp "\$CLEANED" ${haplotype_id}.cleaned.fasta
  """
}