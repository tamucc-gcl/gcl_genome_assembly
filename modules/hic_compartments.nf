process HIC_COMPARTMENTS {
  tag "${meta.id}:${stage}:${resolution}"
  label 'hic_compartments'

  publishDir "${params.outdir}/compartments_plot", mode: params.publish_dir_mode

  input:
    tuple val(meta), val(stage), path(mcool)
    val resolution
    val min_contig_bp
    val max_contigs
    path(compartments_script)

  output:
    tuple val(meta), val(stage),
          path("${meta.id}_${stage}.comp_${resolution}.pc1.bedGraph"),
          emit: pc1_track

    tuple val(meta), val(stage),
          path("${meta.id}_${stage}.comp_${resolution}.pc1.genomewide.png"),
          emit: pc1_png

    tuple val(meta), val(stage),
          path("${meta.id}_${stage}.comp_${resolution}.eigvals.tsv"),
          emit: eigvals

  stub:
  """
  touch ${meta.id}_${stage}.comp_${resolution}.pc1.bedGraph
  touch ${meta.id}_${stage}.comp_${resolution}.pc1.genomewide.png
  touch ${meta.id}_${stage}.comp_${resolution}.eigvals.tsv
  """

  script:
    def out_prefix = "${meta.id}_${stage}.comp_${resolution}"
    def assembly_label = "${meta.id} (${stage})"
  """
  set -euo pipefail

  python ${compartments_script} \\
    --mcool ${mcool} \\
    --resolution ${resolution} \\
    --assembly_id "${assembly_label}" \\
    --out_prefix "${out_prefix}" \\
    --min_contig_bp ${min_contig_bp} \\
    --max_contigs ${max_contigs}
  """
}
