process PANGENOME_REGIONAL_BUILD {
    tag 'pinned-regional-grouping'
    label 'rgfa_patch_build'
    container null
    conda "${projectDir}/environments/rgfa_split_patch.yaml"
    cpus 2
    memory '4 GB'
    time '1h'
    publishDir "${params.outdir}/pangenome/tools/regional_patch", mode: params.publish_dir_mode,
        pattern: 'patched/*', saveAs: { fn -> fn.endsWith('.json') || fn.endsWith('.tsv') || fn.endsWith('.diff') ? fn.replaceFirst('^patched/', '') : null }
    input:
    val(trigger)
    path(sources)
    path(adapter)
    path(adapters)
    output:
    tuple path('patched/rgfa-split-fixed'), path('patched/build.json'), emit: tool
    path('patched/regional-grouping.diff'), emit: patch
    path('patched/native_regressions.tsv'), emit: tests
    script:
    """
    set -euo pipefail
    python3 ${adapter} build --source ${sources} --out patched
    """
}

process PANGENOME_REGIONAL_ASSIGN_FIXED {
    tag "${taxid}:${policy}:regional-fixed"
    label 'cactus_tools'
    cpus 1
    memory '16 GB'
    time '2h'
    publishDir "${params.outdir}/pangenome/${taxid}/deletion_filter", mode: params.publish_dir_mode,
        pattern: 'assignment_*'
    input:
    tuple val(taxid), val(policy), path(gfa), path(paf), path(config), path(selection),
        path(identities), path(provenance), val(contigs), path(binary), path(build_info)
    path(adapter)
    path(adapters)
    output:
    tuple val(taxid), val(policy), path("assignment_${policy}_regional_1Mb_fixed"), emit: results
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
    python3 ${adapter} assign --prepared prepared --binary ${binary} --build-info ${build_info} \\
        --contigs ${contigs.join(' ')} --out assignment_${policy}_regional_1Mb_fixed
    """
}

process PANGENOME_REGIONAL_PATCH_REPORT {
    tag "${taxid}"
    label 'cactus_tools'
    cpus 1
    memory '4 GB'
    time '30m'
    publishDir "${params.outdir}/pangenome/${taxid}/deletion_filter", mode: params.publish_dir_mode,
        pattern: 'summary/*', saveAs: { it.replaceFirst('^summary/', '') }
    input:
    tuple val(taxid), path(original), path(fixed), path(filters)
    path(adapter)
    path(adapters)
    output:
    tuple val(taxid), path('summary/regional_assignment_report.md'), emit: report
    tuple val(taxid), path('summary/*'), emit: files
    script:
    """
    set -euo pipefail
    python3 ${adapter} report --original ${original} --fixed ${fixed} --filters ${filters} --out summary
    """
}

workflow PANGENOME_REGIONAL_PATCH_TEST {
    take:
    filtered
    native_assignments
    filter_reviews
    sources
    adapter
    adapters
    main:
    PANGENOME_REGIONAL_BUILD(filtered.map { row -> true }.unique(), sources, adapter, adapters)
    PANGENOME_REGIONAL_ASSIGN_FIXED(filtered.combine(PANGENOME_REGIONAL_BUILD.out.tool), adapter, adapters)
    PANGENOME_REGIONAL_PATCH_REPORT(native_assignments.filter { taxid, policy, name, dir -> name == 'regional_1Mb' }
        .map { taxid, policy, name, dir -> tuple(taxid,dir) }.groupTuple()
        .join(PANGENOME_REGIONAL_ASSIGN_FIXED.out.results.map { taxid, policy, dir -> tuple(taxid,dir) }.groupTuple())
        .join(filter_reviews.map { taxid, policy, dir -> tuple(taxid,dir) }.groupTuple()), adapter, adapters)
    emit:
    report = PANGENOME_REGIONAL_PATCH_REPORT.out.report
    files = PANGENOME_REGIONAL_PATCH_REPORT.out.files
}
