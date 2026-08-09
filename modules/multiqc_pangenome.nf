/*
========================================================================================
    PANGENOME MultiQC MODULE  (reporting / workstream C)
========================================================================================
    Repo location: modules/multiqc_pangenome.nf

    Aggregates the pangenome QC inputs MultiQC understands into a single report:
      - odgi graph statistics  (odgi module; *.og.stats.yaml -- whole-genome + per chrom)
      - variant catalog stats  (bcftools module; `bcftools stats` on the filtered VCF)

    Both are dropped in the task dir and picked up by `multiqc .`. Additional MultiQC-native
    inputs (e.g. BUSCO-on-graph once wired) can be added to the collected input list.

    Input : tuple(taxid, [mqc_files])   (odgi YAMLs + bcftools stats, collected per taxid)
    Output: report (html) / data / versions
========================================================================================
*/

process MULTIQC_PANGENOME {
    tag "taxid_${taxid}"
    label 'multiqc'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(mqc_files)

    output:
    tuple val(taxid), path("${taxid}_pangenome_multiqc.html"),      emit: report
    tuple val(taxid), path("${taxid}_pangenome_multiqc_data"),      emit: data, optional: true
    path("versions.tsv"),                                           emit: versions

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"

    multiqc --force \\
        --filename ${taxid}_pangenome_multiqc.html \\
        --title "Pangenome QC: ${taxid}" \\
        --comment "odgi graph statistics and bcftools variant-catalog stats for the ${taxid} pangenome" \\
        .

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tmultiqc\\t%s\\n' "${task.process}" "\$(multiqc --version 2>&1 | awk '{print \$NF}')"
    } > versions.tsv
    """

    stub:
    """
    touch ${taxid}_pangenome_multiqc.html
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
