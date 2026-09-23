/*
========================================================================================
    CHIMERIC SCAFFOLD SUBWORKFLOW
========================================================================================
    Repo location: workflows/chimera.nf

    Detection, confirmation and repair of scaffolds that span two reference chromosomes,
    extracted VERBATIM from main.nf.

    WHY IT MOVED
    ------------
    main.nf hit Groovy's compiled-unit limit:

        String too long. The given string is 65676 Unicode code units long, but only a
        maximum of 65535 is allowed.

    The chimera arm is the natural unit to lift out: it is self-contained, it is the newest
    addition, and it accounts for ~9.3 kB of main.nf. Nothing in the logic changed -- the
    bodies below are the same statements in the same order, so behaviour and task hashes are
    unaffected except where a channel had to be renamed to a `take:` parameter.

    WHAT IT DOES
    ------------
      CHIMERA_JOINS      AGP joins + reference PAF -> which joins separate two chromosomes
      CHIMERA_EVIDENCE   Hi-C cross-contact, telomeres, N-gaps -> confirmation per cut
      BREAK_CHIMERAS     applies the concordance gate and cuts, N joins -> N+1 pieces

    `ch_pre_finalize` is emitted rather than assigned in place, because FINALIZE_ASSEMBLY in
    main.nf consumes it. Its default -- harmonized assemblies mixed with the short-read-only
    branch -- is built here so the assign-then-override pattern stays intact; that .mix()
    carries assemblies with no harmonization name map and losing it silently dropped
    short-read-only samples from finalization for several runs.
========================================================================================
*/

include { BREAK_CHIMERAS }   from '../modules/break_chimeras.nf'
include { CHIMERA_JOINS }    from '../modules/chimera_joins.nf'
include { CHIMERA_EVIDENCE } from '../modules/chimera_evidence.nf'

workflow CHIMERA {

    take:
    ch_harmonized                  //  HARMONIZE_SCAFFOLDS.out.assemblies    tuple(meta, fa, nm)
    ch_shortread_finished          //  short-read-only assemblies            tuple(meta, fa)
    ch_round1_agp                  //  SCAFFOLD_HIC_ROUND1.out.agp           tuple(meta, agp)
    ch_round2_agp                  //  round-2 AGP or empty                  tuple(meta, agp)
    ch_ref_pafs_by_id              //  HARMONIZE_SCAFFOLDS.out.ref_pafs_by_id
    ch_chimera_candidates          //  HARMONIZE_SCAFFOLDS.out.chimera_candidates
    ch_ref_name_map                //  HARMONIZE_SCAFFOLDS.out.ref_name_map
    ch_contig_pairs_in             //  FILTER_HIC_BAM.out.pairs (contig stage)
    ch_telo_by_taxid               //  per-taxid telomere motif              tuple(taxid, motif)

    main:
    ch_versions = Channel.empty()
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
    if( params.chimera_detect != false ) {
        ch_agp_script = Channel.fromPath("${projectDir}/py_scripts/agp_joins.py",
                                        checkIfExists: true)
        ch_cj_script  = Channel.fromPath("${projectDir}/py_scripts/chimera_joins.py",
                                        checkIfExists: true)

        // round 2 is conditional, so this join tolerates its absence: without a round-2 AGP
        // the round-1 objects are final, which agp_joins.py handles natively.
        ch_r1_agp = ch_round1_agp.map { meta, agp -> tuple(meta.id, meta.taxid.toString(), agp) }
        ch_r2_agp = ch_round2_agp.map { meta, agp -> tuple(meta.id, agp) }

        CHIMERA_JOINS(
            ch_r1_agp
                .join( ch_r2_agp, remainder: true )
                .join( ch_ref_pafs_by_id, remainder: true )
                .filter { id, taxid, r1, r2, paf -> r1 != null }
                .map { id, taxid, r1, r2, paf ->
                    tuple(taxid, id, r1, r2 ?: file("${projectDir}/assets/NO_ROUND2", checkIfExists: true), paf ?: file("${projectDir}/assets/NO_PAF", checkIfExists: true)) }
                .combine(ch_chimera_candidates, by: 0)
                .combine(ch_ref_name_map, by: 0)
                .map { taxid, id, r1, r2, paf, cand, rnm ->
                    tuple(taxid, id, r1, r2, paf, cand, rnm) },
            ch_agp_script.first(),
            ch_cj_script.first() )
        ch_versions = ch_versions.mix(CHIMERA_JOINS.out.versions)
        ch_chimeric_joins = CHIMERA_JOINS.out.called

        // ---- independent confirmation of each called join ----------------------------
        // Telomere and N-gap evidence, plus a Hi-C cross-contact profile built by
        // TRANSLATING the published contig-space pairs into scaffold coordinates -- no
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

            // the contig-stage pairs: deduplicated UU pairs in contig coordinates
            ch_contig_pairs = ch_contig_pairs_in
                .filter { meta, stage, pairs_gz -> stage == 'contig' }
                .map    { meta, stage, pairs_gz -> tuple(meta.id, pairs_gz) }

            CHIMERA_EVIDENCE(
                CHIMERA_JOINS.out.called
                    .map { taxid, id, called -> tuple(id, taxid, called) }
                    .join( ch_harmonized
                               .map { meta, fa, nm -> tuple(meta.id, fa) } )
                    .join( ch_round1_agp.map { meta, agp -> tuple(meta.id, agp) } )
                    .join( ch_round2_agp.map { meta, agp -> tuple(meta.id, agp) },
                           remainder: true )
                    .join( ch_contig_pairs, remainder: true )
                    .filter { id, taxid, called, fa, r1, r2, pairs -> called != null && fa != null }
                    .map { id, taxid, called, fa, r1, r2, pairs ->
                        tuple(taxid, id, fa, called, r1,
                              r2 ?: file("${projectDir}/assets/NO_ROUND2", checkIfExists: true), pairs ?: file("${projectDir}/assets/NO_PAIRS", checkIfExists: true)) }
                    // the motif is per species, so it attaches by key
                    .combine( ch_telo_by_taxid, by: 0 )
                    .map { taxid, id, fa, called, r1, r2, pairs, motif ->
                        tuple(taxid, id, fa, called, r1, r2, pairs, motif) },
                ch_ce_hic_script.first(),
                ch_ce_script.first() )
            ch_versions = ch_versions.mix(CHIMERA_EVIDENCE.out.versions)
        }
    }

    if( params.chimera_break && params.chimera_break.toString() != 'false' ) {
        ch_break_script = Channel.fromPath("${projectDir}/py_scripts/break_chimeras.py",
                                          checkIfExists: true)
        // 'auto' uses the joins CHIMERA_JOINS just called; a path uses that file, so an
        // edited copy is how you choose which scaffolds to cut. Either way the input is the
        // CALLED JOINS table, not the candidates: it carries cut_bp from the AGP and one row
        // per chimeric join, so a scaffold with N of them yields N+1 pieces.
        ch_cand = ( params.chimera_break.toString() == 'auto' )
            ? ch_chimeric_joins.map { taxid, id, f -> f }.collectFile(name: 'called.tsv',
                                                                      keepHeader: true,
                                                                      skip: 1)
            : Channel.fromPath(params.chimera_break.toString(), checkIfExists: true).first()

        BREAK_CHIMERAS(
            ch_harmonized
                .map { meta, fa, nm -> tuple(meta, fa, nm) }
                .combine( ch_cand ),
            ch_break_script.first() )
        ch_versions = ch_versions.mix(BREAK_CHIMERAS.out.versions)
        ch_broken       = BREAK_CHIMERAS.out.assemblies
        ch_pre_finalize = BREAK_CHIMERAS.out.assemblies
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

    if (params.chimera_detect != false) {
        ch_chimera_candidates_rpt = ch_chimera_candidates
            .map { taxid, f -> f }
            .first()
        // one table per assembly -> one table for the report. keepHeader because every file
        // carries the same header row.
        ch_chimera_joins_rpt = ch_chimeric_joins
            .map { taxid, id, f -> f }
            .collectFile(name: 'all_chimeric_joins.tsv', keepHeader: true, skip: 1)
            .ifEmpty(file('NO_CHIMERA_JOINS'))
    }
    if (params.chimera_detect != false && params.chimera_evidence != false) {
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


    emit:
    pre_finalize        = ch_pre_finalize          // tuple(meta, fasta, name_map)
    broken              = ch_broken                // ONLY the cut assemblies, for QC_PHASE
    called              = ch_chimeric_joins        // tuple(taxid, id, chimeric_joins.tsv)
    candidates_for_rpt  = ch_chimera_candidates_rpt
    joins_for_rpt       = ch_chimera_joins_rpt
    evidence_for_rpt    = ch_chimera_evidence_rpt
    figures_for_rpt     = ch_chimera_figures_rpt
    versions            = ch_versions
}
