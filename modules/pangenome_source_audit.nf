/* Read-only accounting against cached construction products; never rebuilds a graph. */
process PANGENOME_SOURCE_AUDIT {
    tag "${taxid}"
    label 'cactus_tools'
    cpus 1
    memory '4 GB'
    time '2h'
    publishDir "${params.outdir}/pangenome/${taxid}/representation", mode: params.publish_dir_mode,
        pattern: 'source_audit/*', saveAs: { it.replaceFirst('^source_audit/', '') }
    input:
    tuple val(taxid), path(identities), path(paths), path(indexes), val(specs), path(mappings), val(evidence)
    path(adapter)
    output:
    tuple val(taxid), path('source_audit/source_report.md'), emit: report
    tuple val(taxid), path('source_audit/*'), emit: files
    script:
    def config = groovy.json.JsonOutput.toJson([assemblies: specs, evidence: evidence])
    """
    set -euo pipefail
    cat > source_inputs.json <<'SOURCE_INPUTS'
    ${config}
    SOURCE_INPUTS
    python3 ${adapter} --input source_inputs.json --identities ${identities} \\
        --paths ${paths} --out source_audit
    """
}
