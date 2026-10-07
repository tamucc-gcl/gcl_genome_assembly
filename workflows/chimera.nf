include { RECOVER_OLDER_JOINS } from '../modules/recover_older_joins.nf'
include { CHIMERA_COORDINATE_SUMMARY } from '../modules/chimera_coordinate_summary.nf'
/* Chimera assessment and optional cuts on the last-scaffold FASTA, before finishing.
 * The round1 parameter carries the LAST-round AGP; no two-round chain is used.
 * Pairs were mapped to that AGP's input. Older joins across corrections are unresolved.
 * Name assignments are applied only after finishing by FINALIZE_ASSEMBLY.
 */

include { BREAK_CHIMERAS }   from '../modules/break_chimeras.nf'
include { CHIMERA_JOINS }    from '../modules/chimera_joins.nf'
include { CHIMERA_EVIDENCE } from '../modules/chimera_evidence.nf'
include { CHIMERA_SEQUENCE_CONTEXT } from '../modules/chimera_sequence_context.nf'
include { CHIMERA_REVIEW; CHIMERA_REVIEW_INDEX } from '../modules/chimera_review.nf'
include { HARMONIZE_SPECIES as CHIMERA_REASSIGN_SPECIES } from '../modules/harmonize_species.nf'
include { harmonizerArgs } from './harmonize_scaffolds.nf'

workflow CHIMERA {

    take:
    ch_harmonized                  //  HARMONIZE_SCAFFOLDS.out.assemblies    tuple(meta, fa, nm)
    ch_shortread_finished          //  short-read-only assemblies            tuple(meta, fa)
    ch_round1_agp                  // LAST-round AGP, paired with its own input mappings
    ch_ref_pafs_by_id              //  HARMONIZE_SCAFFOLDS.out.ref_pafs_by_id
    ch_chimera_candidates          //  HARMONIZE_SCAFFOLDS.out.chimera_candidates
    ch_ref_name_map                //  HARMONIZE_SCAFFOLDS.out.ref_name_map
    ch_reference_id                // same reference before and after correction
    ch_peer_quality                // measured harmonization voter/passenger status
    ch_contig_pairs_in             // pairs mapped to the LAST-round scaffolding input
    ch_hic_evidence_inputs         // source FASTA and read-set/library identity
    ch_native_graphs               // current native primary graphs, keyed by sample
    ch_telo_by_taxid               //  per-taxid telomere motif              tuple(taxid, motif)

    ch_original_scaffolds
    ch_original_agp
    ch_hifi_reads
    capabilities

    main:
    ch_versions = Channel.empty()
    def retainedBams = params.chimera_hifi_bam_manifest ?
        new groovy.json.JsonSlurper().parseText(file(params.chimera_hifi_bam_manifest, checkIfExists:true).text) : [:]
    def detect_on = capabilities.scaffold && capabilities.harmonize &&
                    params.harmonize_scaffold_names && params.chimera_detect != false
    def evidence_on = detect_on && params.chimera_evidence != false
    if (params.chimera_break?.toString() == 'auto')
        error 'Automated chimera cutting is deferred. Review the evidence and supply an edited chimera_review.tsv.'
    // The broken assemblies on their own, separate from pre_finalize (which mixes them with
    // the untouched and short-read-only ones). main.nf stages THIS as the 'chimera_broken'
    // QC checkpoint, so contiguity before and after a cut is comparable in the QC table.
    ch_broken   = Channel.empty()

    // OFF by default. Run 1 writes <species>.chimera_candidates.tsv and cuts nothing; you
    // review it and supply selected joins. Votes alone never authorize automatic cuts.
    //
    // When off the process is not instantiated and the original channel flows straight
    // through -- no pass-through task and no cache churn.
    // Assigned BEFORE the branch, then overridden inside it. A variable assigned only
    // inside if/else blocks in a workflow body is not visible after them -- "No such
    // variable: ch_pre_finalize". The same assign-then-override pattern is used for
    // ch_hap_priv and ch_hap_cov in workflows/pangenome.nf.
    //
    // The FULL statement: the .mix() carries the short-read-only branch, which has no
    // harmonization name map. An earlier anchor matched only the first line of this and
    // inserted the chimera block between the two, orphaning the .mix() -- which Groovy
    // accepts as a no-op expression, so nothing failed until an output of the never-invoked
    // BREAK_CHIMERAS was read further down.
    ch_pre_finalize = ch_harmonized
        .mix( ch_shortread_finished.map { meta, fa -> tuple(meta, fa, file("${projectDir}/assets/NO_HARMONIZE", checkIfExists: true)) } )

    // ---- which joins are chimeric, and exactly where ---------------------------------
    // Detection runs whenever harmonization did, independently of chimera_break: the
    // candidates and the called joins are the evidence a break is justified by, and they are
    // worth having even on a run that cuts nothing.
    ch_chimeric_joins = Channel.empty()
    ch_evidence_calls = Channel.empty()
    ch_older_summaries = Channel.empty()
    if( detect_on ) {
        ch_agp_script = Channel.fromPath("${projectDir}/py_scripts/agp_joins.py",
                                        checkIfExists: true)
        ch_cj_script  = Channel.fromPath("${projectDir}/py_scripts/chimera_joins.py",
                                        checkIfExists: true)

        // Optional PAFs are a lookup value, not an outer join with unknown empty arity.
        ch_paf_lookup = ch_ref_pafs_by_id.toList()
            .map { records -> records.collectEntries { id, pf -> [(id.toString()): pf] } }
        CHIMERA_JOINS(
            ch_round1_agp.map { meta, agp -> tuple(meta.taxid.toString(), meta.id, agp) }
                .combine(ch_paf_lookup)
                .map { taxid, id, agp, pafs ->
                    tuple(taxid, id, agp,
                        file("${projectDir}/assets/NO_ROUND2", checkIfExists: true),
                        pafs[id.toString()] ?: file("${projectDir}/assets/NO_PAF", checkIfExists: true)) }
                .combine(ch_harmonized.map { meta, fa, nm -> tuple(meta.taxid.toString(), meta.id, fa) }, by: [0, 1])
                .combine(ch_chimera_candidates, by: 0)
                .combine(ch_ref_name_map, by: 0)
                .map { taxid, id, r1, r2, paf, fa, cand, rnm ->
                    tuple(taxid, id, r1, r2, paf, cand, rnm, fa) },
            ch_agp_script.first(),
            ch_cj_script.first(),
            file("${projectDir}/py_scripts/chimera_coordinate_guard.py", checkIfExists: true),
            file("${projectDir}/py_scripts/chimera_intervals.py", checkIfExists: true),
            file("${projectDir}/py_scripts/chimera_schema.py", checkIfExists: true) )
        ch_versions = ch_versions.mix(CHIMERA_JOINS.out.versions)
        ch_chimeric_joins = CHIMERA_JOINS.out.called
        ch_evidence_calls = ch_chimeric_joins
        if (params.run_scaffold_round2 && params.chimera_recover_older) {
            ch_original = ch_original_scaffolds.join(ch_original_agp)
                .map { meta, fa, agp -> tuple(meta.id.toString(), fa, agp) }
            ch_recovery = ch_chimeric_joins
                .map { taxid, id, called -> tuple(id.toString(), taxid, called) }
                .join(ch_original)
                .join(ch_harmonized.map { meta, fa, nm -> tuple(meta.id.toString(), fa) })
                .combine(ch_paf_lookup)
                .map { id, taxid, called, oldfa, oldagp, current, pafs ->
                    tuple(taxid, id, oldfa, oldagp, current, called,
                          pafs[id.toString()] ?: file("${projectDir}/assets/NO_PAF", checkIfExists: true)) }
                .combine(ch_chimera_candidates, by: 0)
                .combine(ch_ref_name_map, by: 0)
            RECOVER_OLDER_JOINS(ch_recovery,
                file("${projectDir}/py_scripts/recover_older_joins.py", checkIfExists: true),
                file("${projectDir}/py_scripts/chimera_schema.py", checkIfExists: true))
            ch_chimeric_joins = RECOVER_OLDER_JOINS.out.called
            ch_evidence_calls = ch_chimeric_joins
            ch_older_summaries = RECOVER_OLDER_JOINS.out.summary.map { taxid, id, report ->
                new groovy.json.JsonSlurper().parseText(report.text)
            }
            ch_versions = ch_versions.mix(RECOVER_OLDER_JOINS.out.versions)
        }

        // Direct local peer comparison is diagnostic only, independent of cutting.
        if (params.chimera_sequence_context) {
            ch_context_quality = ch_peer_quality.toList().map { files ->
                files.collectMany { quality ->
                    def lines = quality.readLines().findAll { it && !it.startsWith('#') }
                    def header = lines[0].split('\t') as List
                    lines.drop(1).collect { line ->
                        def fields = line.split('\t', -1)
                        [(fields[header.indexOf('id')]): fields[header.indexOf('role')] == 'voter' &&
                            !fields[header.indexOf('role_reasons')].toLowerCase().contains('forced')] }
                }.collectEntries() }
            ch_context_cohort = ch_harmonized.toList().map { records -> [records:records] }
                .combine(ch_context_quality)
                .map { cohort, quality -> [records: cohort.records.collect { m, fa, nm ->
                    def labels = [:]
                    if (nm.name != 'NO_HARMONIZE') {
                        def lines = nm.readLines().findAll { it && !it.startsWith('#') }
                        def header = lines[0].split('\t') as List
                        lines.drop(1).each { line ->
                            def fields = line.split('\t', -1)
                            def name = fields[header.indexOf('new_name')]
                            if (header.contains('chromosome_member') && fields[header.indexOf('chromosome_member')]=='yes' && name.startsWith('chr') && !name.contains('+'))
                                labels[fields[header.indexOf('old_name')]] = name.split('_')[0]
                        }
                    }
                    tuple(m + [auto_evidence: quality[m.id.toString()] == true, chromosome_labels:labels], fa, nm) }] }
            ch_context_reads = ch_hifi_reads.toList()
                .map { records -> records.collectEntries { m, fq -> [(m.sample.toString()): fq] } }
            ch_context_graphs = ch_native_graphs.toList().map { records ->
                records.collectEntries { m, paths -> [(m.sample.toString()): paths instanceof List ? paths : [paths]] } }
            ch_context_assays = ch_hic_evidence_inputs
                .map { m, fa, libraries -> tuple(m.id.toString(), fa, libraries) }
                .join(ch_round1_agp.map { m, agp -> tuple(m.id.toString(), agp) })
                .join(ch_contig_pairs_in.map { m, stage, pairs -> tuple(m.id.toString(), pairs) })
                .toList()
                .map { records -> records.collectEntries { id, source, libraries, agp, pairs ->
                    [(id): [source:source, libraries:libraries, agp:agp, pairs:pairs]] } }
            ch_context_in = ch_chimeric_joins
                .map { taxid, id, calls -> tuple(id.toString(), calls) }
                .join(ch_harmonized.map { m, fa, nm -> tuple(m.id.toString(), m, fa) })
                .combine(ch_context_cohort)
                .combine(ch_context_reads)
                .combine(ch_context_assays)
                .combine(ch_context_graphs)
                .map { id, calls, m, fa, cohort, reads, assays, graphs ->
                    def peers = cohort.records.findAll { pm, pf, pn -> pm.taxid.toString() == m.taxid.toString() && pm.id != m.id }
                    def assay = assays[id] ?: [:]
                    def absent = file("${projectDir}/assets/NO_PAIRS", checkIfExists: true)
                    def retained = retainedBams[id]
                    def suffix = m.id.toString().startsWith(m.sample.toString()+'_') ? m.id.toString().substring(m.sample.toString().length()+1) : 'primary'
                    def nativeGraph = graphs[m.sample.toString()]?.find { it.name == m.sample.toString()+'.'+suffix+'.p_ctg.gfa' }
                    tuple(m.taxid.toString(), m, fa, calls,
                          peers.collect { it[0] }, peers ? peers.collect { it[1] } : [file("${projectDir}/assets/NO_PAF", checkIfExists: true)],
                          !retained && params.chimera_hifi_context && reads[m.sample.toString()] ? reads[m.sample.toString()] : absent,
                          assay.agp ?: absent, assay.pairs ?: absent, assay.source ?: absent, assay.libraries ?: absent,
                          retained ? [file(retained.bam,checkIfExists:true),file(retained.index,checkIfExists:true)] : [absent],
                          retained ? file(retained.provenance,checkIfExists:true) : absent, nativeGraph ?: absent)
                }
                .combine(ch_telo_by_taxid, by: 0)
                .map { taxid, m, fa, calls, pm, pf, reads, agp, pairs, source, libraries, bams, provenance, graph, motif ->
                    tuple(m, fa, calls, pm, pf, reads, motif ?: 'CCCTAA', agp, pairs, source, libraries, bams, provenance, graph) }
            CHIMERA_SEQUENCE_CONTEXT(ch_context_in,
                Channel.value(['chimera_sequence_context.py', 'chimera_controls.py', 'chimera_graph_evidence.py', 'chimera_blocks.py', 'chimera_tracks.py'].collect {
                    file("${projectDir}/py_scripts/${it}", checkIfExists: true) }))
            ch_versions = ch_versions.mix(CHIMERA_SEQUENCE_CONTEXT.out.versions)
            CHIMERA_REVIEW(
                ch_chimeric_joins.map { taxid, id, calls -> tuple(id.toString(), calls) }
                    .join(CHIMERA_SEQUENCE_CONTEXT.out.context.map { m, context -> tuple(m.id.toString(), m, context) })
                    .map { id, calls, m, context -> tuple(m, calls, context) },
                Channel.value(['chimera_review.py', 'chimera_markdown.py'].collect { file("${projectDir}/py_scripts/${it}", checkIfExists: true) }))
            CHIMERA_REVIEW_INDEX(CHIMERA_REVIEW.out.packets.map { m, packet -> packet }.toList(),
                file("${projectDir}/py_scripts/chimera_review_index.py", checkIfExists: true))
        }

        // ---- supplementary evidence for each called join ----------------------------
        // Telomere and N-gap evidence, plus a Hi-C cross-contact profile built by
        // TRANSLATING last-round input pairs into current scaffold coordinates -- no
        // re-alignment, and no dependence on a contact map that does not exist yet at this
        // point in the DAG.
        //
        // Runs on a NON-BREAKING run by design: the evidence is what justifies a cut, so it
        // must exist before one is made. After a break it would look for an interstitial
        // array on a scaffold that no longer exists.
        if( params.chimera_evidence != false ) {
            ch_ce_hic_script = Channel.fromPath("${projectDir}/py_scripts/chimera_hic_pairs.py",
                                               checkIfExists: true)
            ch_ce_script     = Channel.fromPath("${projectDir}/py_scripts/chimera_evidence.py",
                                               checkIfExists: true)

            // Deduplicated pairs aligned to the last scaffolding input.
            ch_contig_pairs = ch_contig_pairs_in
                .map    { meta, stage, pairs_gz -> tuple(meta.id, pairs_gz) }

            CHIMERA_EVIDENCE(
                ch_evidence_calls
                    .map { taxid, id, called -> tuple(id.toString(), taxid, called) }
                    .join( ch_harmonized
                               .map { meta, fa, nm -> tuple(meta.id.toString(), fa) } )
                    .join( ch_round1_agp.map { meta, agp -> tuple(meta.id.toString(), agp) } )
                    .combine(ch_contig_pairs.toList()
                        .map { records -> records.collectEntries { id, pf -> [(id.toString()): pf] } })
                    .map { id, taxid, called, fa, agp, pairs ->
                        tuple(taxid, id, fa, called, agp,
                              file("${projectDir}/assets/NO_ROUND2", checkIfExists: true),
                              pairs[id.toString()] ?: file("${projectDir}/assets/NO_PAIRS", checkIfExists: true)) }
                    // the motif is per species, so it attaches by key
                    .combine( ch_telo_by_taxid, by: 0 )
                    .map { taxid, id, fa, called, r1, r2, pairs, motif ->
                        tuple(taxid, id, fa, called, r1, r2, pairs, motif) },
                ch_ce_hic_script.first(),
                ch_ce_script.first() )
            ch_versions = ch_versions.mix(CHIMERA_EVIDENCE.out.versions)
        }
    }

    if( capabilities.scaffold && capabilities.harmonize && params.harmonize_scaffold_names && params.chimera_break && params.chimera_break.toString() != 'false' ) {
        ch_break_script = Channel.value(['break_chimeras.py', 'chimera_actions.py']
            .collect { file("${projectDir}/py_scripts/${it}", checkIfExists: true) })
        // Review file is a staged input: changing it invalidates cutting/downstream tasks only.
        ch_known_assemblies = ch_harmonized.map { m, fa, nm -> m.id.toString() }.toList().map { ids -> [known:ids] }
        ch_review_file = Channel.fromPath(params.chimera_break.toString(), checkIfExists: true).first()
            .combine(ch_known_assemblies)
            .map { review, catalog ->
                def lines = review.readLines().findAll { it && !it.startsWith('#') }
                def header = lines[0].split('\t') as List
                if (!header.contains('selected') || !header.contains('assembly')) error 'Review file requires selected and assembly columns'
                lines.drop(1).each { line ->
                    def values = line.split('\t',-1)
                    def selection = values[header.indexOf('selected')]
                    if (!(selection in ['YES','NO'])) error 'Every review row requires selected=YES or NO'
                    if (selection=='YES' && !(values[header.indexOf('assembly')] in catalog.known))
                        error 'Selected review row names an unknown assembly: '+values[header.indexOf('assembly')]
                }
                review }
        ch_break_in = ch_harmonized.combine(ch_review_file)
        ch_break_in.branch { meta, fa, nm, cand ->
            def lines = cand.readLines().findAll { it && !it.startsWith('#') }
            def header = lines[0].split('\t') as List
            def selected = lines.drop(1).any { line ->
                def row = line.split('\t',-1)
                row[header.indexOf('selected')]=='YES' && row[header.indexOf('assembly')]==meta.id.toString() }
            cut: nm.name != 'NO_HARMONIZE' && selected
            passthrough: true
        }.set { ch_break_routes }
        BREAK_CHIMERAS(ch_break_routes.cut, ch_break_script.first())
        ch_versions = ch_versions.mix(BREAK_CHIMERAS.out.versions)
        ch_broken = BREAK_CHIMERAS.out.assemblies.map { m,fa,nm -> tuple(m.id.toString(),m,fa,nm) }
            .join(BREAK_CHIMERAS.out.verification.map { m,verification -> tuple(m.id.toString(),verification) })
            .filter { id,m,fa,nm,verification -> new groovy.json.JsonSlurper().parseText(verification.text).selected_actions.size() > 0 }
            .map { id,m,fa,nm,verification -> tuple(m,fa,nm) }
        ch_pre_finalize = BREAK_CHIMERAS.out.assemblies
            .mix(ch_break_routes.passthrough.map { meta, fa, nm, cand -> tuple(meta, fa, nm) })
            .mix(ch_shortread_finished.map { meta, fa -> tuple(meta, fa, file("${projectDir}/assets/NO_HARMONIZE", checkIfExists: true)) })
        // Reassign from actual corrected sequences, never from parent side labels.
        ch_reassign_lr = ch_pre_finalize.filter { m, fa, nm -> m.assembler == 'hifiasm' }
        ch_changed_taxids = BREAK_CHIMERAS.out.verification
            .filter { m, verification -> new groovy.json.JsonSlurper().parseText(verification.text).selected_actions.size() > 0 }
            .map { m, verification -> m.taxid.toString() }.toList().map { ids -> [taxids:ids] }
        ch_reassign_species = ch_reassign_lr.map { m, fa, nm -> tuple(m.taxid.toString(), m.id.toString(), fa) }
            .groupTuple()
            .filter { taxid, ids, fastas -> ids.size() >= 2 }
            .map { taxid, ids, fastas ->
                def order = (0..<ids.size()).toList().sort { ids[it] }
                tuple(taxid, order.collect { ids[it] }, order.collect { fastas[it] }) }
            .combine(ch_changed_taxids)
            .filter { taxid, ids, fastas, changed -> taxid in changed.taxids }
            .map { taxid, ids, fastas, changed -> tuple(taxid, ids, fastas) }
            .join(ch_reference_id.map { taxid, id -> tuple(taxid.toString(), id) })
        CHIMERA_REASSIGN_SPECIES(ch_reassign_species,
            file("${projectDir}/py_scripts/harmonize_names.py", checkIfExists: true), harmonizerArgs())
        ch_versions = ch_versions.mix(CHIMERA_REASSIGN_SPECIES.out.versions)
        ch_reassigned_maps = CHIMERA_REASSIGN_SPECIES.out.name_maps
            .flatMap { taxid, maps ->
                (maps instanceof List ? maps : [maps]).collect { nm ->
                    tuple(nm.name.replaceFirst(/\.harmonized_name_map\.tsv$/, ''), nm) } }
        ch_pre_finalize = ch_reassign_lr.map { m, fa, nm -> tuple(m.id.toString(), m, fa, nm) }
            .join(ch_reassigned_maps, remainder: true)
            .filter { row -> row[1] != null }
            .map { id, m, fa, oldMap, newMap -> tuple(m, fa, newMap ?: oldMap) }
            .mix(ch_shortread_finished.map { m, fa -> tuple(m, fa, file("${projectDir}/assets/NO_HARMONIZE", checkIfExists: true)) })
    }

    ch_name_map_files = ch_harmonized
        .map { meta, fa, nm -> nm }

    // ---- chimera tables for the report ---------------------------------------------
    // Defaulted to sentinels BEFORE any branch: assigned only inside one, they would not be
    // visible at the REPORTING call, which is the "No such variable" failure this workflow
    // body has hit repeatedly.
    ch_chimera_candidates_rpt = Channel.value(file('NO_CHIMERA_CANDIDATES'))
    ch_chimera_joins_rpt      = Channel.value(file('NO_CHIMERA_JOINS'))
    ch_chimera_evidence_rpt   = Channel.value(file('NO_CHIMERA_EVIDENCE'))
    ch_chimera_figures_rpt    = Channel.value(file('NO_CHIMERA_FIGURES'))

    if (detect_on) {
        ch_chimera_candidates_rpt = ch_chimera_candidates
            .map { taxid, f -> f }
            .first()
            .ifEmpty(file('NO_CHIMERA_CANDIDATES'))
        // one table per assembly -> one table for the report. keepHeader because every file
        // carries the same header row.
        ch_chimera_joins_rpt = ch_chimeric_joins
            .map { taxid, id, f -> f }
            .collectFile(name: 'all_chimeric_joins.tsv', keepHeader: true, skip: 1)
            .ifEmpty(file('NO_CHIMERA_JOINS'))
    }
    if (evidence_on) {
        // each evidence file is ONE cut in metric/value long form, so they go as a list and
        // the R side widens and stacks them
        ch_chimera_evidence_rpt = CHIMERA_EVIDENCE.out.evidence
            .map { taxid, id, files -> files }
            .flatten()
            .collect()
            .ifEmpty([file('NO_CHIMERA_EVIDENCE')])
        // the figures travel as their own channel: the R side cannot derive them from the
        // evidence paths, because those stage as bare filenames with no sibling .png present
        ch_chimera_figures_rpt = CHIMERA_EVIDENCE.out.figures
            .map { taxid, id, files -> files }
            .flatten()
            .collect()
            .ifEmpty([file('NO_CHIMERA_FIGURES')])
    }


    CHIMERA_COORDINATE_SUMMARY(
        ch_harmonized.map { meta, fa, nm ->
            [id: meta.id, hic: meta.hic, harmonized: nm.name != 'NO_HARMONIZE']
        }.mix(ch_shortread_finished.map { meta, fa ->
            [id: meta.id, hic: meta.hic, harmonized: false]
        }).toList(),
        ch_chimeric_joins.map { taxid, id, f -> id.toString() }.toList(), ch_older_summaries.toList())

    emit:
    coordinate_report = CHIMERA_COORDINATE_SUMMARY.out.report
    pre_finalize        = ch_pre_finalize          // tuple(meta, fasta, name_map)
    broken              = ch_broken                // ONLY the cut assemblies, for QC_PHASE
    called              = ch_chimeric_joins        // tuple(taxid, id, chimeric_joins.tsv)
    candidates_for_rpt  = ch_chimera_candidates_rpt
    joins_for_rpt       = ch_chimera_joins_rpt
    evidence_for_rpt    = ch_chimera_evidence_rpt
    figures_for_rpt     = ch_chimera_figures_rpt
    versions            = ch_versions
}
