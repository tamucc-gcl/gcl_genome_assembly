process COLLECT_SOFTWARE_VERSIONS {
    tag "software_versions"
    conda 'conda-forge::python=3.11'
    cpus 1
    memory '2 GB'
    time '2h'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode, pattern: 'provenance/software_versions.tsv', saveAs: { 'software_versions.tsv' }
    publishDir "${params.outdir}/provenance", mode: params.publish_dir_mode, pattern: 'provenance/*', saveAs: { it.replaceFirst('^provenance/','') }

    input:
    path version_files, stageAs:'versions/version??/*'
    val sources
    val context
    path helper
    path catalogs, stageAs:'catalog/source??/*'

    output:
    path 'provenance/software_versions.tsv', emit: versions
    path 'provenance/*', emit: files

    script:
    """
    cat > sources.json <<'SOURCES'
${groovy.json.JsonOutput.toJson(sources)}
SOURCES
    cat > context.json <<'CONTEXT'
${groovy.json.JsonOutput.toJson(context)}
CONTEXT
    python3 ${helper} --specs sources.json --staged versions --context context.json --out provenance
    """
}
