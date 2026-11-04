process ANNOTATION_AUDIT {
    tag "${meta.sample}"
    cpus 1
    memory '1 GB'
    time '15m'
    conda "${projectDir}/environments/annotation_audit.yaml"
    errorStrategy 'terminate'
    publishDir "${params.outdir}/reports/organelle_annotations", mode: params.publish_dir_mode
    input:
    tuple val(meta),path(genbank)
    path helper
    output:
    tuple val(meta),path("${genbank.simpleName}.counts.tsv"),emit: counts
    tuple val(meta),path("${genbank.simpleName}.loci.tsv"),emit: loci
    script:
    """
    python3 ${helper} --genbank '${genbank}' --sample '${meta.sample}' --prefix '${genbank.simpleName}'
    """
}
