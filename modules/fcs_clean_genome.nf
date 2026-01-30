process FCS_CLEAN_GENOME {
  tag "${haplotype_id}"
  label 'fcs' 
  
  publishDir "${params.outdir}/${stage}/decontam", mode: params.publish_dir_mode

  input:
    tuple val(haplotype_id), path(assembly_fa), path(action_report), val(stage)

  output:
    tuple val(haplotype_id), path("${haplotype_id}.decontaminated.fasta"), emit: decontaminated_fasta
    tuple val(haplotype_id), path("${haplotype_id}.contaminants.fasta"),   emit: contaminants_fasta
    
  script:
  """
  set -euo pipefail

    /app/bin/gx clean-genome \\
    --action-report ${action_report} \\
    --input ${assembly_fa} \\
    --output ${haplotype_id}.decontaminated.fasta \\
    --contam-fasta-out ${haplotype_id}.contaminants.fasta
  """
}