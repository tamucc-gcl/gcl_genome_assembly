/* Production tool build is independent of removable replay/trace workflows. */
process PANGENOME_SEQUENCE_TOOLS {
    tag 'guarded-regional-tools'
    label 'rgfa_patch_build'
    container null
    conda "${projectDir}/environments/rgfa_split_patch.yaml"
    cpus 2
    memory '4 GB'
    time '1h'
    publishDir "${params.outdir}/pangenome/tools/sequence_policy", mode: params.publish_dir_mode,
        pattern: 'sequence_tool_provenance/*', saveAs: { it.replaceFirst('^sequence_tool_provenance/', '') }
    input:
    val(trigger)
    path(source)
    path(adapter)
    path(helpers)
    output:
    path('sequence_tools'), emit: tools
    path('sequence_tool_provenance/build.json'), emit: provenance
    script:
    """
    set -euo pipefail
    python3 ${adapter} build --source ${source} --out sequence_tools
    mkdir sequence_tool_provenance
    cp sequence_tools/build.json sequence_tool_provenance/
    """
    stub:
    """
    mkdir -p sequence_tools sequence_tool_provenance
    printf '{"status":"STUB"}\\n' > sequence_tools/build.json
    cp sequence_tools/build.json sequence_tool_provenance/
    """
}
