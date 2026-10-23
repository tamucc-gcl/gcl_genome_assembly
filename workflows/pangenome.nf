include { PANGENOME_CONSTRUCTION } from './pangenome_construction.nf'
include { PANGENOME_ANALYSIS } from './pangenome_analysis.nf'
include { PANGENOME_VARIANT_AUDIT } from '../modules/pangenome_variant_audit.nf'
include { PANGENOME_BUILD_QC } from '../modules/pangenome_build_qc.nf'
include { PANGENOME_SOURCE_AUDIT } from '../modules/pangenome_source_audit.nf'
include { PANGENOME_ASSIGNMENT_PREPARE; PANGENOME_ASSIGNMENT_REPLAY; PANGENOME_ASSIGNMENT_REPORT } from '../modules/pangenome_assignment_replay.nf'
include { PANGENOME_DELETION_DIAGNOSTIC } from '../modules/pangenome_deletion_replay.nf'

/* Species selection is owned by the always-on assembly eligibility audit. */
workflow PANGENOME {
    take:
    ch_eligibility_manifest
    main:
    ch_versions = Channel.empty()
    ch_report = Channel.empty()
    ch_manifest = Channel.empty()
    if (params.run_pangenome) {
        if (params.pangenome_regional_chunk_replay && !params.pangenome_regional_gate_trace)
            error 'Regional chunk comparison requires pangenome_regional_gate_trace'
        if (params.pangenome_regional_gate_trace && !params.pangenome_regional_patch_replay)
            error 'Regional gate tracing requires pangenome_regional_patch_replay'
        if (params.pangenome_regional_patch_replay && (!params.pangenome_deletion_filter_replay || params.pangenome_validate_only))
            error 'Regional patch replay requires the deletion-filter diagnostic and a completed graph (not validation-only)'
        if (params.pangenome_gref != true || params.pangenome_gref_source != 'clip')
            error 'The default graph contract requires pangenome_gref=true and pangenome_gref_source=clip'
        if (params.pangenome_grefl != null)
            error 'Allele-clustered GREF is outside the standard catalog contract; leave pangenome_grefl unset'
        ch_cohorts = ch_eligibility_manifest.flatMap { manifest ->
            def data = new groovy.json.JsonSlurper().parseText(manifest.text)
            data.cohorts.findAll { it.ready }.collect { cohort ->
                tuple(cohort.taxid, cohort.reference_name,
                    cohort.members.collect { it.graph_name },
                    cohort.members.collect { file(it.fasta, checkIfExists: true) },
                    cohort.reference_contigs,
                    cohort + [min_individuals: data.settings.min_individuals])
            }
        }
        PANGENOME_CONSTRUCTION(ch_cohorts)
        ch_versions = PANGENOME_CONSTRUCTION.out.versions
        ch_report = PANGENOME_CONSTRUCTION.out.report
        if (!params.pangenome_validate_only) {
            PANGENOME_ANALYSIS(PANGENOME_CONSTRUCTION.out.clip_stats_input,
                PANGENOME_CONSTRUCTION.out.clip_sharing_input, ch_eligibility_manifest)
            ch_report = ch_report.mix(PANGENOME_ANALYSIS.out.report)
            PANGENOME_VARIANT_AUDIT(PANGENOME_CONSTRUCTION.out.variant_input,
                file("${projectDir}/py_scripts/pangenome_variant_audit.py", checkIfExists: true))
            ch_report = ch_report.mix(PANGENOME_VARIANT_AUDIT.out.report)
            PANGENOME_BUILD_QC(PANGENOME_CONSTRUCTION.out.support_qc
                .join(PANGENOME_ANALYSIS.out.stats)
                .join(PANGENOME_CONSTRUCTION.out.identities)
                .map { taxid, bundle, full, versions, stats, identities ->
                    tuple(taxid, bundle, full, stats, identities, versions)
                },
                file("${projectDir}/py_scripts/pangenome_build_qc.py", checkIfExists: true))
            ch_report = ch_report.mix(PANGENOME_BUILD_QC.out.report)
            ch_source_inputs = ch_eligibility_manifest.flatMap { manifest ->
                def data = new groovy.json.JsonSlurper().parseText(manifest.text)
                data.cohorts.findAll { it.ready }.collect { cohort ->
                    def members = cohort.members.collect { member -> data.assemblies.find { it.id == member.id } }
                    def indexes = members.collect { file(it.fai, checkIfExists: true) }
                    def specs = members.withIndex().collect { member, i ->
                        [id: member.id, fai: indexes[i].name,
                         chromosomes: member.chromosome_scale ? member.reference_contigs : []]
                    }
                    tuple(cohort.taxid, indexes, specs)
                }
            }
            ch_source_mappings = PANGENOME_CONSTRUCTION.out.source_evidence
            PANGENOME_SOURCE_AUDIT(PANGENOME_BUILD_QC.out.files
                .map { taxid, files -> tuple(taxid, files.find { it.name == 'full_clip_sequence_paths.tsv' }) }
                .join(PANGENOME_CONSTRUCTION.out.identities).join(ch_source_inputs).join(ch_source_mappings)
                .map { taxid, paths, identities, indexes, specs, mappings, evidence ->
                    tuple(taxid, identities, paths, indexes, specs, mappings, evidence)
                }, file("${projectDir}/py_scripts/pangenome_source_audit.py", checkIfExists: true))
            ch_report = ch_report.mix(PANGENOME_SOURCE_AUDIT.out.report)
            if (params.pangenome_assignment_replay) {
                if (params.pangenome_cactus_extra)
                    error 'Assignment replay currently requires the pinned default assignment configuration (no pangenome_cactus_extra)'
                def assignment_adapters = ['pangenome_assignment_replay.py','pangenome_source_audit.py']
                    .collect { file("${projectDir}/py_scripts/${it}", checkIfExists: true) }
                PANGENOME_ASSIGNMENT_PREPARE(PANGENOME_CONSTRUCTION.out.assignment_input
                    .join(PANGENOME_CONSTRUCTION.out.identities)
                    .join(PANGENOME_SOURCE_AUDIT.out.files.map { taxid, files ->
                        tuple(taxid, files.find { it.name == 'source_sequences.tsv' })
                    }), assignment_adapters)
                PANGENOME_ASSIGNMENT_REPLAY(PANGENOME_ASSIGNMENT_PREPARE.out.prepared.flatMap { taxid, graph, paf, config, selection, identities, provenance, contigs ->
                    ['baseline','whole_relaxed','regional_1Mb','regional_3Mb'].collect { name ->
                        tuple(taxid,graph,paf,config,selection,identities,provenance,contigs,name)
                    }
                }, assignment_adapters)
                PANGENOME_ASSIGNMENT_REPORT(PANGENOME_ASSIGNMENT_REPLAY.out.results.groupTuple()
                    .map { taxid, names, directories -> tuple(taxid,directories) }
                    .join(PANGENOME_ASSIGNMENT_PREPARE.out.coverage), assignment_adapters)
                ch_report = ch_report.mix(PANGENOME_ASSIGNMENT_REPORT.out.report)
            }
            if (params.pangenome_deletion_filter_replay) {
                if (params.pangenome_cactus_extra)
                    error 'Deletion replay requires the pinned default configuration (no pangenome_cactus_extra)'
                def deletion_adapters = ['pangenome_deletion_replay.py','pangenome_assignment_replay.py','pangenome_source_audit.py']
                    .collect { file("${projectDir}/py_scripts/${it}", checkIfExists: true) }
                PANGENOME_DELETION_DIAGNOSTIC(PANGENOME_CONSTRUCTION.out.deletion_input
                    .join(PANGENOME_CONSTRUCTION.out.identities)
                    .join(PANGENOME_SOURCE_AUDIT.out.files.map { taxid, files ->
                        tuple(taxid, files.find { it.name == 'source_sequences.tsv' })
                    }), deletion_adapters,
                    file("${projectDir}/vendor/cactus-gfa-tools-regional", checkIfExists: true),
                    file("${projectDir}/py_scripts/pangenome_regional_patch.py", checkIfExists: true))
                ch_report = ch_report.mix(PANGENOME_DELETION_DIAGNOSTIC.out.report)
            }
            // The legacy construction `cactus --version` row reports Toil, not Cactus.
            // Keep it in the build provenance but use the separately captured versions here.
            ch_versions = PANGENOME_ANALYSIS.out.versions.mix(PANGENOME_VARIANT_AUDIT.out.versions,
                                                           PANGENOME_BUILD_QC.out.versions)
        }
        ch_manifest = PANGENOME_CONSTRUCTION.out.manifest
    }
    emit:
    report = ch_report
    manifest = ch_manifest
    versions = ch_versions
}
