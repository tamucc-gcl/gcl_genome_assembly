/* Optional cached-input experiment. No construction task or assembly is changed. */
include { PANGENOME_REGIONAL_PATCH_TEST } from './pangenome_regional_patch.nf'
process PANGENOME_DELETION_PREPARE {
    tag "${taxid}"
    label 'cactus_tools'
    cpus 6
    memory '32 GB'
    time '4h'
    publishDir "${params.outdir}/pangenome/${taxid}/deletion_filter", mode: params.publish_dir_mode,
        pattern: 'preparation_review'
    input:
    tuple val(taxid), path(gfa), path(paf), path(unfiltered), path(filter_log), val(contigs), path(identities), path(sources)
    path(adapters)
    output:
    tuple val(taxid), path('deletion_prepared'), val(contigs), emit: prepared
    tuple val(taxid), path('preparation_review'), emit: provenance
    script:
    """
    set -euo pipefail
    python3 pangenome_deletion_replay.py prepare --gfa ${gfa} --paf ${paf} \\
        --unfiltered ${unfiltered} --filter-log ${filter_log} --identities ${identities} \\
        --sources ${sources} --threads ${task.cpus} --out deletion_prepared
    mkdir preparation_review
    cp deletion_prepared/{selection.tsv,identities.tsv,effective_config.xml,preparation.json,original_filter.log} preparation_review/
    """
}

process PANGENOME_DELETION_FILTER {
    tag "${taxid}:${policy}"
    label 'cactus_tools'
    cpus 6
    memory '32 GB'
    time '4h'
    publishDir "${params.outdir}/pangenome/${taxid}/deletion_filter", mode: params.publish_dir_mode,
        pattern: 'review_*'
    input:
    tuple val(taxid), path(prepared), val(contigs), val(policy)
    path(adapters)
    output:
    tuple val(taxid), val(policy), path("policy_${policy}/cactus.sv.gfa"), path("policy_${policy}/filtered.paf"),
        path("policy_${policy}/effective_config.xml"), path("policy_${policy}/selection.tsv"),
        path("policy_${policy}/identities.tsv"), path("policy_${policy}/preparation.json"), val(contigs), emit: prepared
    tuple val(taxid), val(policy), path("review_${policy}"), emit: reviews
    script:
    """
    set -euo pipefail
    python3 pangenome_deletion_replay.py filter --prepared ${prepared} --policy ${policy} \\
        --threads ${task.cpus} --out policy_${policy}
    mkdir review_${policy}
    for source in policy_${policy}/*.json policy_${policy}/*.tsv policy_${policy}/*.tsv.gz policy_${policy}/*.log; do
        [[ ! -f "\$source" ]] || cp "\$source" review_${policy}/
    done
    """
}

process PANGENOME_DELETION_ASSIGN {
    tag "${taxid}:${policy}:${case_name}"
    label 'cactus_tools'
    cpus 1
    memory '16 GB'
    time '2h'
    publishDir "${params.outdir}/pangenome/${taxid}/deletion_filter", mode: params.publish_dir_mode,
        pattern: 'assignment_*'
    input:
    tuple val(taxid), val(policy), path(gfa), path(paf), path(config), path(selection),
        path(identities), path(provenance), val(contigs), val(case_name)
    path(adapters)
    output:
    tuple val(taxid), val(policy), val(case_name), path("assignment_${policy}_${case_name}"), emit: results
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
    python3 pangenome_assignment_replay.py replay --prepared prepared --case '${case_name}' \\
        --contigs ${contigs.join(' ')} --out assignment_${policy}_${case_name}
    printf '{"policy":"${policy}","case":"${case_name}"}\\n' > assignment_${policy}_${case_name}/case.json
    """
}

process PANGENOME_DELETION_REPORT {
    tag "${taxid}"
    label 'cactus_tools'
    cpus 1
    memory '4 GB'
    time '30m'
    publishDir "${params.outdir}/pangenome/${taxid}/deletion_filter", mode: params.publish_dir_mode,
        pattern: 'summary/*', saveAs: { it.replaceFirst('^summary/', '') }
    input:
    tuple val(taxid), path(filters), path(assignments)
    path(adapters)
    output:
    tuple val(taxid), path('summary/deletion_filter_report.md'), emit: report
    tuple val(taxid), path('summary/*'), emit: files
    script:
    """
    set -euo pipefail
    python3 pangenome_deletion_replay.py report --filters ${filters} \\
        --assignments ${assignments} --out summary
    """
}

workflow PANGENOME_DELETION_DIAGNOSTIC {
    take:
    ch_inputs
    adapters
    patch_sources
    patch_adapter
    main:
    PANGENOME_DELETION_PREPARE(ch_inputs, adapters)
    PANGENOME_DELETION_FILTER(PANGENOME_DELETION_PREPARE.out.prepared.flatMap { taxid, prepared, contigs ->
        ['native','chromosome_isolated','disabled'].collect { policy -> tuple(taxid,prepared,contigs,policy) }
    }, adapters)
    PANGENOME_DELETION_ASSIGN(PANGENOME_DELETION_FILTER.out.prepared.flatMap { taxid, policy, graph, paf, config, selection, identities, provenance, contigs ->
        ['baseline','regional_1Mb'].collect { name -> tuple(taxid,policy,graph,paf,config,selection,identities,provenance,contigs,name) }
    }, adapters)
    PANGENOME_DELETION_REPORT(PANGENOME_DELETION_FILTER.out.reviews.groupTuple()
        .map { taxid, policies, dirs -> tuple(taxid,dirs) }
        .join(PANGENOME_DELETION_ASSIGN.out.results.groupTuple()
            .map { taxid, policies, cases, dirs -> tuple(taxid,dirs) }), adapters)
    ch_report = PANGENOME_DELETION_REPORT.out.report
    if (params.pangenome_regional_patch_replay) {
        PANGENOME_REGIONAL_PATCH_TEST(PANGENOME_DELETION_FILTER.out.prepared,
            PANGENOME_DELETION_ASSIGN.out.results, PANGENOME_DELETION_FILTER.out.reviews,
            patch_sources, patch_adapter, adapters)
        ch_report = ch_report.mix(PANGENOME_REGIONAL_PATCH_TEST.out.report)
    }
    emit:
    report = ch_report
    files = PANGENOME_DELETION_REPORT.out.files
}
