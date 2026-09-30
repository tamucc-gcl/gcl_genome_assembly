include { PANGENOME_CONSTRUCTION } from './pangenome_construction.nf'
include { PANGENOME_ANALYSIS } from './pangenome_analysis.nf'

/* Species selection is owned by the always-on assembly eligibility audit. */
workflow PANGENOME {
    take:
    ch_eligibility_manifest
    main:
    ch_versions = Channel.empty()
    ch_report = Channel.empty()
    ch_manifest = Channel.empty()
    if (params.run_pangenome) {
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
        if (!params.pangenome_validate_only) {
            PANGENOME_ANALYSIS(PANGENOME_CONSTRUCTION.out.clip_stats_input)
            ch_versions = ch_versions.mix(PANGENOME_ANALYSIS.out.versions)
        }
        ch_report = PANGENOME_CONSTRUCTION.out.report
        ch_manifest = PANGENOME_CONSTRUCTION.out.manifest
    }
    emit:
    report = ch_report
    manifest = ch_manifest
    versions = ch_versions
}
