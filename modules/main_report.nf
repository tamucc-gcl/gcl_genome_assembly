process MAIN_REPORT {
    tag 'main_report'
    label 'main_reporting'
    errorStrategy 'terminate'
    publishDir "${params.outdir}", mode: params.publish_dir_mode
    input:
    path eligibility
    val expected
    path name_maps,stageAs:'names/*'
    path paintings,stageAs:'painting/*'
    path qc_metrics,stageAs:'qc/*'
    path decisions,stageAs:'decisions/*'
    path verifications,stageAs:'verification/*'
    val artifacts
    val context
    path versions,stageAs:'versions/*'
    path renderer
    output:
    path 'assembly_report.md',emit: report
    path 'reports/main',emit: assets
    script:
    def spec=[eligibility:eligibility.toString(),expected:expected,artifacts:artifacts,run:context,
        name_maps:name_maps.collect { it.toString() },paintings:paintings.collect { it.toString() },
        qc:qc_metrics.toString(),decisions:decisions.toString(),
        verifications:verifications.collect { it.toString() },versions:versions.toString()]
    """
    set -euo pipefail
    export TZ="\${TZ:-UTC}"
    cat > report_spec.json <<'REPORT_SPEC'
${groovy.json.JsonOutput.toJson(spec)}
REPORT_SPEC
    Rscript ${renderer} report_spec.json
    """
}
