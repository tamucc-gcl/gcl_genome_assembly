/*
========================================================================================
    ORGANELLE_ANNOTATION — annotate assembled organelles (selector subworkflow)
========================================================================================
    Repo location: workflows/organelle_annotation.nf

    Sibling to SHORTREAD_ORGANELLE; consumes its .out.assembly (+ .out.stats for status) and
    routes each organelle to the right annotator:

      animal_mt / fungus_mt  -> MITOS2 (ANNOTATE_MITO), refseq89m / refseq89f, --linear if not circular
      embplant_pt            -> Chloë  (ANNOTATE_PLASTID)
      embplant_mt            -> NO annotator (no maintained CLI tool); emit a note pointing to
                                GeSeq (web) for the report

    Genetic code for the mito branch comes from ch_gcode (geneticCodeFor); plastid code (11)
    is handled inside Chloë. Circular-vs-linear status is joined in from ch_stats.

    Owns its DB/setup (DOWNLOAD_MITOS_DB, SETUP_CHLOE), gated on whether the relevant organelle
    types are actually present — same self-contained pattern as SHORTREAD_ORGANELLE.

    take:
      ch_assembly  tuple(meta, org_type, fasta)     SHORTREAD_ORGANELLE.out.assembly
      ch_stats     tuple(meta, org_type, stats_tsv) SHORTREAD_ORGANELLE.out.stats
      ch_gcode     tuple(taxid, genetic_code)       mito genetic code per taxon
    emit:
      mito_annotation     tuple(meta, org_type, bed)
      plastid_annotation  tuple(meta, org_type, gff)
      plant_mito_note     tuple(meta, org_type, note)   GeSeq pointer for the report
      versions
========================================================================================
*/

include { DOWNLOAD_MITOS_DB } from '../modules/download_mitos_db.nf'
include { SETUP_CHLOE }       from '../modules/setup_chloe.nf'
include { ANNOTATE_MITO }     from '../modules/annotate_mito.nf'
include { ANNOTATE_PLASTID }  from '../modules/annotate_plastid.nf'
include { mitosRefseqFor }    from '../functions/taxonomy.nf'

// Plant mitochondrion: no maintained CLI annotator. Emit a note (picked up by the report)
// telling the user to annotate manually with GeSeq.
process PLANT_MITO_NOTE {
    tag "${meta.sample}:${org_type}"
    label 'tiny'
    publishDir "${params.outdir}/organelle/annotation/${org_type}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(org_type)

    output:
    tuple val(meta), val(org_type), path("${meta.sample}.${org_type}.annotation_note.txt"), emit: note

    script:
    """
    cat > "${meta.sample}.${org_type}.annotation_note.txt" <<'EOF'
Plant mitochondrion assembled but NOT automatically annotated.
There is no maintained command-line annotator for plant mitogenomes; annotate manually with
GeSeq (web-only): https://chlorobox.mpimp-golm.mpg.de/geseq.html
EOF
    """

    stub:
    """
    echo "stub note" > ${meta.sample}.${org_type}.annotation_note.txt
    """
}

workflow ORGANELLE_ANNOTATION {
    take:
    ch_assembly
    ch_stats
    ch_gcode

    main:
    // Attach status (circular/linear) from stats onto each assembly, keyed by (sample, org_type).
    ch_status = ch_stats.map { meta, org, tsv ->
        def cols = tsv.text.trim().readLines().last().split('\t')
        tuple( [meta.sample, org], (cols.size() > 2 ? cols[2] : 'unknown') )
    }

    ch_asm = ch_assembly
        .map { meta, org, fa -> tuple([meta.sample, org], meta, org, fa) }
        .join( ch_status )
        .map { key, meta, org, fa, status -> tuple(meta, org, status, fa) }

    ch_asm
        .branch { meta, org, status, fa ->
            mito:    org.endsWith('_mt') && org != 'embplant_mt'   // animal_mt, fungus_mt
            plastid: org == 'embplant_pt'
            plantmt: org == 'embplant_mt'
        }
        .set { br }

    // ── metazoan / fungal mito → MITOS2 ──────────────────────────────────────
    ch_mitos_db = DOWNLOAD_MITOS_DB(
                      br.mito.map { meta, org, status, fa -> params.mitos_downloads }.unique(),
                      params.mitos_refseq_sets,
                      params.mitos_force_download
                  ).db.first()

    br.mito
        .map { meta, org, status, fa -> tuple(meta.taxid?.toString(), meta, org, status, fa) }
        .combine( ch_gcode, by: 0 )
        .map { taxid, meta, org, status, fa, gcode ->
               tuple(meta, org, status, gcode, (params.mitos_refseq ?: mitosRefseqFor(org)), fa) }
        .set { ch_mito_in }
    ANNOTATE_MITO( ch_mito_in, ch_mitos_db )

    // ── plant plastid → Chloë ────────────────────────────────────────────────
    ch_chloe = SETUP_CHLOE(
                   br.plastid.map { meta, org, status, fa -> params.chloe_downloads }.unique(),
                   params.chloe_force_setup
               ).dir.first()
    ANNOTATE_PLASTID( br.plastid.map { meta, org, status, fa -> tuple(meta, org, fa) }, ch_chloe )

    // ── plant mito → no annotation, GeSeq note for the report ────────────────
    PLANT_MITO_NOTE( br.plantmt.map { meta, org, status, fa -> tuple(meta, org) } )

    emit:
    mito_annotation    = ANNOTATE_MITO.out.bed
    plastid_annotation = ANNOTATE_PLASTID.out.gff
    plant_mito_note    = PLANT_MITO_NOTE.out.note
    versions           = ANNOTATE_MITO.out.versions.mix( ANNOTATE_PLASTID.out.versions )
}
