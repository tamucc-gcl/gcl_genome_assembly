/*
========================================================================================
    RIPARIAN (RIBBON) SYNTENY PLOT MODULE
========================================================================================
    Generates a riparian / ribbon synteny plot from a minimap2 PAF alignment
    between two haplotype assemblies using gggenomes.

    Input:
    - haplotype_id1: Reference haplotype identifier (top row)
    - assembly1:     Reference assembly FASTA
    - haplotype_id2: Query haplotype identifier (bottom row)
    - assembly2:     Query assembly FASTA
    - paf:           Gzipped PAF alignment from PAIRWISE_ALIGNMENT

    Output:
    - riparian: Riparian plot PNG

    Dependencies: R with gggenomes (CRAN), ggplot2, dplyr, argparse, RColorBrewer
========================================================================================
*/

process RIPARIAN_PLOT {
    tag "${haplotype_id1}_vs_${haplotype_id2}"
    label 'pairwise_alignment'

    publishDir "${params.outdir}/pairwise_alignments",
        mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id1), path(assembly1), val(haplotype_id2), path(assembly2), path(paf)

    output:
    tuple val(haplotype_id1), val(haplotype_id2), path("${haplotype_id1}_vs_${haplotype_id2}_riparian.png"), emit: riparian

    script:
    def min_aln   = params.riparian_min_aln_bp   ?: 50000
    def min_seq   = params.riparian_min_seq_bp    ?: 1000000
    def alpha     = params.riparian_alpha         ?: 0.45
    def width     = params.riparian_width         ?: 14
    def height    = params.riparian_height        ?: 6
    def out_prefix = "${haplotype_id1}_vs_${haplotype_id2}"
    """
    set -euo pipefail

    # Generate FAI indices for both assemblies
    samtools faidx ${assembly1}
    samtools faidx ${assembly2}

    # Generate riparian plot
    Rscript ${projectDir}/r_scripts/riparian_paf.R \\
        --paf       "${paf}" \\
        --ref_fai   "${assembly1}.fai" \\
        --query_fai "${assembly2}.fai" \\
        --ref       "${haplotype_id1}" \\
        --query     "${haplotype_id2}" \\
        --output    "${out_prefix}_riparian.png" \\
        --width     ${width} \\
        --height    ${height} \\
        --min_aln   ${min_aln} \\
        --min_seq   ${min_seq} \\
        --alpha     ${alpha}
    """

    stub:
    def out_prefix = "${haplotype_id1}_vs_${haplotype_id2}"
    """
    touch "${out_prefix}_riparian.png"
    """
}