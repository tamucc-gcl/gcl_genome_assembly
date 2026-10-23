include { PANGENOME_REGIONAL_CHUNK_TEST } from './pangenome_regional_chunks.nf'

process PANGENOME_REGIONAL_TRACE_BUILD {
    tag 'regional-observation'
    label 'rgfa_patch_build'
    container null
    conda "${projectDir}/environments/rgfa_split_patch.yaml"
    cpus 2
    memory '4 GB'
    time '1h'
    publishDir "${params.outdir}/pangenome/tools/regional_trace", mode: params.publish_dir_mode,
        pattern: 'trace_provenance/*', saveAs: { it.replaceFirst('^trace_provenance/', '') }
    input:
    val(trigger)
    path(sources)
    path(adapter)
    path(helpers)
    output:
    path('trace_tools'), emit: tools
    path('trace_provenance/*'), emit: provenance
    script:
    """
    set -euo pipefail
    python3 ${adapter} build --source ${sources} --out trace_tools
    mkdir trace_provenance
    cp trace_tools/trace_build.json trace_tools/*.diff trace_tools/gate_regressions.tsv trace_provenance/
    """
}

process PANGENOME_REGIONAL_TRACE {
    tag "${taxid}:${policy}:gates"
    label 'cactus_tools'
    cpus 1
    memory '16 GB'
    time '4h'
    publishDir "${params.outdir}/pangenome/${taxid}/deletion_filter", mode: params.publish_dir_mode,
        pattern: 'trace_*'
    input:
    tuple val(taxid), val(policy), path(gfa), path(paf), path(config), path(selection),
        path(identities), path(provenance), val(contigs), path(binary), path(build_info), path(tools)
    path(adapter)
    path(helpers)
    output:
    tuple val(taxid), val(policy), path("trace_${policy}"), emit: results
    script:
    """
    set -euo pipefail
    mkdir prepared
    ln -s "\$PWD/${gfa}" prepared/cactus.sv.gfa
    ln -s "\$PWD/${paf}" prepared/filtered.paf
    ln -s "\$PWD/${config}" prepared/effective_config.xml
    ln -s "\$PWD/${selection}" prepared/selection.tsv
    ln -s "\$PWD/${identities}" prepared/identities.tsv
    ln -s "\$PWD/${provenance}" prepared/preparation.json
    python3 ${adapter} run --prepared prepared --tools ${tools} --fixed-binary ${binary} --fixed-info ${build_info} \\
        --contigs ${contigs.join(' ')} --out trace_${policy}
    """
}

process PANGENOME_REGIONAL_TRACE_REPORT {
    tag "${taxid}"
    label 'cactus_tools'
    cpus 1
    memory '4 GB'
    time '30m'
    publishDir "${params.outdir}/pangenome/${taxid}/deletion_filter", mode: params.publish_dir_mode,
        pattern: 'gate_summary/*', saveAs: { it.replaceFirst('^gate_summary/', '') }
    input:
    tuple val(taxid), path(traces)
    path(adapter)
    path(helpers)
    output:
    tuple val(taxid), path('gate_summary/regional_gate_report.md'), emit: report
    tuple val(taxid), path('gate_summary/*'), emit: files
    script:
    """
    set -euo pipefail
    python3 ${adapter} report --inputs ${traces} --out gate_summary
    """
}

workflow PANGENOME_REGIONAL_TRACE_TEST {
    take:
    filtered
    fixed_tool
    sources
    adapter
    helpers
    main:
    PANGENOME_REGIONAL_TRACE_BUILD(filtered.map { row -> true }.unique(), sources, adapter, helpers)
    PANGENOME_REGIONAL_TRACE(filtered.combine(fixed_tool).combine(PANGENOME_REGIONAL_TRACE_BUILD.out.tools), adapter, helpers)
    PANGENOME_REGIONAL_TRACE_REPORT(PANGENOME_REGIONAL_TRACE.out.results
        .map { taxid, policy, dir -> tuple(taxid,dir) }.groupTuple(), adapter, helpers)
    ch_reports = PANGENOME_REGIONAL_TRACE_REPORT.out.report
    if (params.pangenome_regional_chunk_replay) {
        PANGENOME_REGIONAL_CHUNK_TEST(filtered, fixed_tool, PANGENOME_REGIONAL_TRACE_BUILD.out.tools,
            PANGENOME_REGIONAL_TRACE.out.results,
            file("${projectDir}/py_scripts/pangenome_regional_chunks.py", checkIfExists: true), helpers + [adapter])
        ch_reports = ch_reports.mix(PANGENOME_REGIONAL_CHUNK_TEST.out.report)
    }
    emit:
    report = ch_reports
}
