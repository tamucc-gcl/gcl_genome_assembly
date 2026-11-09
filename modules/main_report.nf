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
    // A singleton staged path is a Path, which Groovy can iterate by path
    // component. Normalize it before collect (verification/NO_PAIRS is ONE file).
    def fileNames = { value ->
        (value instanceof Collection ? value : [value]).collect { it.toString() }
    }
    def spec=[eligibility:eligibility.toString(),expected:expected,artifacts:artifacts,run:context,
        name_maps:fileNames(name_maps),paintings:fileNames(paintings),
        qc:qc_metrics.toString(),decisions:decisions.toString(),
        verifications:fileNames(verifications),versions:fileNames(versions)]
    """
    set -euo pipefail
    export TZ="\${TZ:-UTC}"
    cat > report_spec.json <<'REPORT_SPEC'
${groovy.json.JsonOutput.toJson(spec)}
REPORT_SPEC
    Rscript ${renderer} report_spec.json
    """
}
