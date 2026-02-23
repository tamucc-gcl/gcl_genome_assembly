process HIC_COMPARTMENTS {
  tag "${haplotype_id}:${stage}:${resolution}"
  label 'hic_compartments'

  publishDir "${params.outdir}/compartments_plot", mode: params.publish_dir_mode

  input:
    tuple val(haplotype_id), val(stage), path(mcool)
    val resolution
    val min_contig_bp
    val max_contigs

  output:
    tuple val(haplotype_id), val(stage),
          path("${haplotype_id}_${stage}.comp_${resolution}.pc1.bedGraph"),
          emit: pc1_track

    tuple val(haplotype_id), val(stage),
          path("${haplotype_id}_${stage}.comp_${resolution}.pc1.genomewide.png"),
          emit: pc1_png

    tuple val(haplotype_id), val(stage),
          path("${haplotype_id}_${stage}.comp_${resolution}.eigvals.tsv"),
          emit: eigvals

  script:
    def out_prefix = "${haplotype_id}_${stage}.comp_${resolution}"
    def assembly_label = "${haplotype_id} (${stage})"
  """
  set -euo pipefail

  ${projectDir}/py_scripts/plot_compartments_pc1_genomewide.py \\
    --mcool ${mcool} \\
    --resolution ${resolution} \\
    --assembly_id "${assembly_label}" \\
    --out_prefix "${out_prefix}" \\
    --min_contig_bp ${min_contig_bp} \\
    --max_contigs ${max_contigs}
  """
}