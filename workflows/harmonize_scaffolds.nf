/*
========================================================================================
    HARMONIZE SCAFFOLDS SUBWORKFLOW
========================================================================================
    Repo location: workflows/harmonize_scaffolds.nf

    Reference-guided scaffold-name harmonization across same-species assemblies, run on
    the pre-FINALIZE channel. Long-read assemblies of a species with >= 2 assemblies are
    harmonized against one in-batch reference; short-read assemblies and single-assembly
    species pass through untouched. Emits tuple(meta, fasta, name_map|NO_HARMONIZE) for
    FINALIZE_ASSEMBLY, which applies the map or falls back to size-rank naming.

    Grouping is by meta.taxid, so mixed-species runs harmonize each species independently
    (e.g. short-read-only samples finalize normally while long-read samples harmonize).

    Toggle: params.harmonize_scaffold_names (default true). When false the whole thing is
    a pass-through with no gather/barrier.

    REFERENCE SELECTION. With params.harmonize_two_pass_selection (default true) every
    plausible candidate is scored by actually aligning against it and building the join
    graph, then the best is used -- see modules/harmonize_two_pass.nf. Set it false to fall
    back to the in-module bash heuristic (count nearest the batch median, N50 tie-break).
    A pinned params.harmonize_reference_ids short-circuits both.

    The harmonizer argument string is built ONCE here and passed to both HARMONIZE_SCORE
    and HARMONIZE_SPECIES, so a candidate is always scored under exactly the settings its
    frame would be built with. Duplicating the argument construction per process would let
    the two drift apart silently.

    Take : ch_assemblies = tuple(meta, fasta)   (the ch_final_assembly mix, pre-FINALIZE)
           resolver       = py_scripts/harmonize_names.py
    Emit : assemblies = tuple(meta, fasta, name_map|NO_HARMONIZE)
           report / versions
========================================================================================
*/

include { HARMONIZE_SPECIES } from '../modules/harmonize_species.nf'
include { HARMONIZE_CANDIDATES; HARMONIZE_SCORE; HARMONIZE_SELECT } from '../modules/harmonize_two_pass.nf'


// Every harmonize_names.py flag shared by scoring and the real run. Built once so the two
// cannot diverge; per-process flags (--outdir, --score-only, --list-candidates) are added
// by the processes themselves.
def harmonizerArgs() {
    def p = params
    def nz = { k, d -> (p[k] != null) ? p[k] : d }          // 0 is falsy in Groovy
    def a = []
    a << "--min-scaffold-bp ${nz('harmonize_min_scaffold_bp', nz('finalize_min_scaffold_bp', 1000000))}"
    a << "--chromosome-set-method ${p.harmonize_chromosome_method ?: 'dropoff'}"
    a << "--dropoff-ratio ${nz('harmonize_dropoff_ratio', 2.0)}"
    a << "--dropoff-min-frac ${nz('harmonize_dropoff_min_frac', 0.5)}"
    a << "--min-aligned-frac ${nz('harmonize_min_aligned_frac', 0.5)}"
    a << "--contained-frac ${nz('harmonize_contained_frac', 0.9)}"
    a << ((p.harmonize_demote_contained == false) ? '--no-demote-contained' : '--demote-contained')
    a << "--secondary-frac ${nz('harmonize_secondary_frac', 0.2)}"
    a << "--min-chrom-frac ${nz('harmonize_min_chrom_frac', 0.1)}"
    a << ((p.harmonize_batch_consensus == false) ? '--no-batch-consensus' : '--batch-consensus')
    a << "--batch-consensus-action ${p.harmonize_batch_consensus_action ?: 'flag'}"
    a << "--voter-min-cut-ratio ${nz('harmonize_voter_min_cut_ratio', 0.0)}"
    a << "--voter-min-genome-frac ${nz('harmonize_voter_min_genome_frac', 0.8)}"
    a << "--voter-min-n50-ratio ${nz('harmonize_voter_min_n50_ratio', 0.2)}"
    a << ((p.harmonize_voter_require_dropoff == false) ? '--no-voter-require-dropoff' : '--voter-require-dropoff')
    a << "--min-voters ${nz('harmonize_min_voters', 3)}"
    a << "--restricted-presence-min-individuals ${nz('harmonize_restricted_presence_min_individuals', 2)}"
    if( p.harmonize_concordance_exclude_reference == true ) a << '--concordance-exclude-reference'
    a << "--member-cover-frac ${nz('harmonize_member_cover_frac', 0.5)}"
    a << "--inflated-aln-tol ${nz('harmonize_inflated_aln_tol', 1.05)}"
    a << "--graph-rule ${p.harmonize_graph_rule ?: 'majority'}"
    a << "--graph-min-fused ${nz('harmonize_graph_min_fused', 2)}"
    return a.join(' ')
}

workflow HARMONIZE_SCAFFOLDS {

    take:
    ch_assemblies    // tuple(meta, fasta)
    resolver         // path to py_scripts/harmonize_names.py

    main:
    ch_versions = Channel.empty()

    // declared here rather than passed from main.nf: main.nf sits ~1.1 kB below Groovy's
    // 65,535-char compiled-unit limit, so lines there are the scarce resource
    selector = file("${projectDir}/py_scripts/select_reference.py", checkIfExists: true)

    if( !params.harmonize_scaffold_names ) {
        // pass-through: every assembly gets the sentinel, no barrier
        ch_out    = ch_assemblies.map { meta, fa -> tuple(meta, fa, file('NO_HARMONIZE')) }
        ch_report = Channel.empty()
        ch_ref_id = Channel.empty()
        ch_ref_scores = Channel.empty()
    }
    else {
        // long-read candidates vs short-read (short-read never harmonized)
        ch_assemblies
            .branch { meta, fa ->
                lr: !meta.shortread
                sr:  meta.shortread
            }
            .set { ch_split }

        // keep an id-keyed copy of the long-read assemblies to rejoin name maps to
        ch_lr_by_id = ch_split.lr.map { meta, fa -> tuple(meta.id, meta, fa) }

        // group long-read assemblies by species; only harmonize groups of >= 2
        ch_species = ch_split.lr
            .map { meta, fa -> tuple(meta.taxid?.toString(), meta.id, fa) }
            .groupTuple()
            .filter { taxid, ids, fas -> ids.size() >= 2 }
            // groupTuple() emits in arrival order, which is not stable across runs. Sort ids
            // and assemblies together with one permutation so ids[i] <-> assemblies[i] is
            // preserved (the module contract) and the task hash stops churning.
            .map { taxid, ids, fas ->
                def perm = (0..<ids.size()).toList().sort { ids[it] }
                tuple(taxid, perm.collect { ids[it] }, perm.collect { fas[it] })
            }

        def hargs = harmonizerArgs()

        // Two-pass only earns its cost when nothing is pinned: a pinned id is already the
        // answer, and scoring candidates to ignore the result would be pure waste.
        def two_pass = (params.harmonize_two_pass_selection != false) \
                       && !(params.harmonize_reference_ids)

        if( two_pass ) {
            HARMONIZE_CANDIDATES( ch_species, resolver, hargs )
            ch_versions = ch_versions.mix( HARMONIZE_CANDIDATES.out.versions )

            // fan out: one scoring task per candidate. readLines() on the emitted table is
            // the only place this workflow turns file CONTENT into channel cardinality, so
            // a parsing slip would silently yield an empty channel rather than an error --
            // hence the ifEmpty guard below and the >= 2 candidate check in the process.
            ch_score_in = ch_species
                .join( HARMONIZE_CANDIDATES.out.candidates )
                .flatMap { taxid, ids, fas, cfile ->
                    cfile.readLines()
                         .findAll { it?.trim() && !it.startsWith('#') && !it.startsWith('candidate\t') }
                         .collect { line -> tuple(taxid, line.tokenize('\t')[0], ids, fas) }
                }
                .ifEmpty { error "HARMONIZE: two-pass selection produced no candidates to score. " +
                                 "Check <taxid>.reference_candidates.tsv, or set " +
                                 "harmonize_two_pass_selection = false." }

            HARMONIZE_SCORE( ch_score_in, resolver, hargs )
            ch_versions = ch_versions.mix( HARMONIZE_SCORE.out.versions.first() )

            HARMONIZE_SELECT( HARMONIZE_SCORE.out.score.groupTuple(), selector )
            ch_versions   = ch_versions.mix( HARMONIZE_SELECT.out.versions )
            ch_ref_scores = HARMONIZE_SELECT.out.scores

            ch_harm_in = ch_species
                .join( HARMONIZE_SELECT.out.reference_id.map { t, f -> tuple(t, f.text.trim()) } )
        }
        else {
            ch_ref_scores = Channel.empty()
            ch_harm_in = ch_species.map { t, i, f -> tuple(t, i, f, 'NO_SELECTION') }
        }

        HARMONIZE_SPECIES( ch_harm_in, resolver, hargs )
        ch_versions = ch_versions.mix( HARMONIZE_SPECIES.out.versions )
        ch_report   = HARMONIZE_SPECIES.out.report

        // re-key emitted name maps by assembly id (filename = <id>.harmonized_name_map.tsv)
        ch_name_map_by_id = HARMONIZE_SPECIES.out.name_maps
            .transpose()
            .map { taxid, nm ->
                def id = nm.name.replaceFirst(/\.harmonized_name_map\.tsv$/, '')
                tuple(id, nm)
            }

        // harmonized long-read assemblies get their map; singletons fall through to sentinel
        ch_lr_out = ch_lr_by_id
            .join( ch_name_map_by_id, remainder: true )
            .map { id, meta, fa, nm -> tuple(meta, fa, nm ?: file('NO_HARMONIZE')) }

        // short-read always sentinel
        ch_sr_out = ch_split.sr.map { meta, fa -> tuple(meta, fa, file('NO_HARMONIZE')) }

        ch_out = ch_lr_out.mix( ch_sr_out )

        ch_ref_id = HARMONIZE_SPECIES.out.reference_id.map { taxid, f -> tuple(taxid, f.text.trim()) }
    }

    emit:
    assemblies     = ch_out      // tuple(meta, fasta, name_map|NO_HARMONIZE)
    reference_id   = ch_ref_id
    reference_scores = ch_ref_scores
    report         = ch_report
    versions       = ch_versions
}
