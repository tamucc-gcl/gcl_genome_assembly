process HIC_TADS {
  tag "${meta.id}:${stage}:${resolution}:w${window_bp}"
  label 'hic_tads'

  input:
    tuple val(meta), val(stage), path(mcool)
    val resolution
    val window_bp
    val min_contig_bp
    val max_contigs

  output:
    tuple val(meta), val(stage),
          path("${meta.id}_${stage}.tads_${resolution}bp_w${window_bp}.insulation.tsv"),
          emit: insulation_tsv

    tuple val(meta), val(stage),
          path("${meta.id}_${stage}.tads_${resolution}bp_w${window_bp}.boundaries.tsv"),
          emit: boundaries_tsv

    tuple val(meta), val(stage),
          path("${meta.id}_${stage}.tads_${resolution}bp_w${window_bp}.tad_book.pdf"),
          emit: tad_book_pdf

  script:
    def out_prefix = "${meta.id}_${stage}.tads_${resolution}bp_w${window_bp}"
    def assembly_label = "${meta.id} (${stage})"
    def balance_flag = params.hic_balance ? "--balance" : ""

  """
  set -euo pipefail

  python ${projectDir}/py_scripts/make_tad_book.py \\
    --mcool ${mcool} \\
    --resolution ${resolution} \\
    --window_bp ${window_bp} \\
    --assembly_id "${assembly_label}" \\
    --out_prefix "${out_prefix}" \\
    --min_contig_bp ${min_contig_bp} \\
    --max_contigs ${max_contigs} \\
    ${balance_flag}

  mv "${out_prefix}.insulation.tsv" "${meta.id}_${stage}.tads_${resolution}bp_w${window_bp}.insulation.tsv"
  mv "${out_prefix}.boundaries.tsv" "${meta.id}_${stage}.tads_${resolution}bp_w${window_bp}.boundaries.tsv"
  mv "${out_prefix}.tad_book.pdf" "${meta.id}_${stage}.tads_${resolution}bp_w${window_bp}.tad_book.pdf"
  """
}
