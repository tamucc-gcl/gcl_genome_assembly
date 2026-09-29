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

workflow CHIMERA {

    take:
    ch_harmonized                  //  HARMONIZE_SCAFFOLDS.out.assemblies    tuple(meta, fa, nm)
    ch_shortread_finished          //  short-read-only assemblies            tuple(meta, fa)
    ch_round1_agp                  // LAST-round AGP, paired with its own input mappings
    ch_ref_pafs_by_id              //  HARMONIZE_SCAFFOLDS.out.ref_pafs_by_id
    ch_chimera_candidates          //  HARMONIZE_SCAFFOLDS.out.chimera_candidates
    ch_ref_name_map                //  HARMONIZE_SCAFFOLDS.out.ref_name_map
    ch_contig_pairs_in             // pairs mapped to the LAST-round scaffolding input
    ch_telo_by_taxid               //  per-taxid telomere motif              tuple(taxid, motif)

    ch_original_scaffolds
    ch_original_agp
    capabilities

    main:
    ch_versions = Channel.empty()
    def detect_on = capabilities.scaffold && capabilities.harmonize &&
                    params.harmonize_scaffold_names && params.chimera_detect != false
    def evidence_on = detect_on && params.chimera_evidence != false
    if (params.chimera_break?.toString() == 'auto' && !detect_on)
        error 'chimera_break=auto requires enabled chimera detection and a harmonization cohort'
    // The broken assemblies on their own, separate from pre_finalize (which mixes them with
    // the untouched and short-read-only ones). main.nf stages THIS as the 'chimera_broken'
    // QC checkpoint, so contiguity before and after a cut is comparable in the QC table.
    ch_broken   = Channel.empty()

    // OFF by default. Run 1 writes <species>.chimera_candidates.tsv and cuts nothing; you
    // review it and pass it back, or set 'auto' to cut what the vote already flagged.
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
        ch_paf_lookup = ch_ref_pafs_by_id.collect(flat: false)
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
            file("${projectDir}/py_scripts/chimera_intervals.py", checkIfExists: true) )
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
                file("${projectDir}/py_scripts/recover_older_joins.py", checkIfExists: true))
            ch_chimeric_joins = RECOVER_OLDER_JOINS.out.called
            ch_evidence_calls = ch_chimeric_joins
            ch_older_summaries = RECOVER_OLDER_JOINS.out.summary.map { taxid, id, report ->
                new groovy.json.JsonSlurper().parseText(report.text)
            }
            ch_versions = ch_versions.mix(RECOVER_OLDER_JOINS.out.versions)
        }

        // ---- independent confirmation of each called join ----------------------------
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
                    .combine(ch_contig_pairs.collect(flat: false)
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
        ch_break_script = Channel.fromPath("${projectDir}/py_scripts/break_chimeras.py",
                                          checkIfExists: true)
        // 'auto' uses the joins CHIMERA_JOINS just called; a path uses that file, so an
        // edited copy is how you choose which scaffolds to cut. Either way the input is the
        // CALLED JOINS table, not the candidates: it carries cut_bp from the AGP and one row
        // per chimeric join, so a scaffold with N of them yields N+1 pieces.
        ch_break_in = Channel.empty()
        if (params.chimera_break.toString() == 'auto') {
            // Left-side pass-through keeps every assembly when no join table was available.
            ch_break_in = ch_harmonized.map { meta, fa, nm -> tuple(meta.id, meta, fa, nm) }
                .join(ch_chimeric_joins.map { taxid, id, f -> tuple(id, f) }, remainder: true)
                .filter { row -> row[1] != null }
                .map { id, meta, fa, nm, cand ->
                    tuple(meta, fa, nm, cand ?: file("${projectDir}/assets/NO_HARMONIZE", checkIfExists: true)) }
        } else {
            ch_break_in = ch_harmonized.combine(
                Channel.fromPath(params.chimera_break.toString(), checkIfExists: true).first())
        }
        ch_break_in.branch { meta, fa, nm, cand ->
            cut: nm.name != 'NO_HARMONIZE' && cand.name != 'NO_HARMONIZE'
            passthrough: true
        }.set { ch_break_routes }
        BREAK_CHIMERAS(ch_break_routes.cut, ch_break_script.first())
        ch_versions = ch_versions.mix(BREAK_CHIMERAS.out.versions)
        ch_broken       = BREAK_CHIMERAS.out.assemblies
        ch_pre_finalize = BREAK_CHIMERAS.out.assemblies
            .mix(ch_break_routes.passthrough.map { meta, fa, nm, cand -> tuple(meta, fa, nm) })
            .mix(ch_shortread_finished.map { meta, fa -> tuple(meta, fa, file("${projectDir}/assets/NO_HARMONIZE", checkIfExists: true)) })
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
        }).collect(flat: false),
        ch_chimeric_joins.map { taxid, id, f -> id.toString() }.collect(), ch_older_summaries.collect(flat: false))

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
