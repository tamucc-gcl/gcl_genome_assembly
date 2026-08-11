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

    Take : ch_assemblies = tuple(meta, fasta)   (the ch_final_assembly mix, pre-FINALIZE)
           resolver       = py_scripts/harmonize_names.py
    Emit : assemblies = tuple(meta, fasta, name_map|NO_HARMONIZE)
           report / versions
========================================================================================
*/

include { HARMONIZE_SPECIES } from '../modules/harmonize_species.nf'

workflow HARMONIZE_SCAFFOLDS {

    take:
    ch_assemblies    // tuple(meta, fasta)
    resolver         // path to py_scripts/harmonize_names.py

    main:
    ch_versions = Channel.empty()

    if( !params.harmonize_scaffold_names ) {
        // pass-through: every assembly gets the sentinel, no barrier
        ch_out    = ch_assemblies.map { meta, fa -> tuple(meta, fa, file('NO_HARMONIZE')) }
        ch_report = Channel.empty()
        ch_ref_id = Channel.empty()
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

        HARMONIZE_SPECIES( ch_species, resolver )
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
    assemblies = ch_out      // tuple(meta, fasta, name_map|NO_HARMONIZE)
    reference_id = ch_ref_id
    report     = ch_report
    versions   = ch_versions
}
