#!/usr/bin/env nextflow
// Standalone test: SHORTREAD_ORGANELLE (assembly) -> ORGANELLE_ANNOTATION (annotation).
// Hands both subworkflows the same channels their parents build; nothing is re-derived.
//   reads sheet : sample,taxid,reads_1,reads_2
//   specs sheet : taxid,type,recursion,kmers,coverage,word_size   (k-list uses ';'; word_size blank => auto)
//   gcode sheet : taxid,genetic_code   (mito annotation code; only taxa with animal_mt/fungus_mt need a row)
nextflow.enable.dsl = 2

include { SHORTREAD_ORGANELLE }  from '../workflows/shortread_organelle.nf'
include { ORGANELLE_ANNOTATION } from '../workflows/organelle_annotation.nf'

params.reads_samplesheet = null
params.organelle_specs   = null
params.gcode_samplesheet = null

workflow {
    if( !params.reads_samplesheet || !params.organelle_specs )
        error "Provide --reads_samplesheet and --organelle_specs (and --gcode_samplesheet to annotate mito)"

    ch_reads = Channel.fromPath( params.reads_samplesheet, checkIfExists: true )
        .splitCsv( header: true )
        .map { row ->
            tuple( [ sample: row.sample, taxid: row.taxid ],
                   file(row.reads_1, checkIfExists: true),
                   file(row.reads_2, checkIfExists: true) )
        }

    ch_organelle = Channel.fromPath( params.organelle_specs, checkIfExists: true )
        .splitCsv( header: true )
        .map { row ->
            tuple( row.taxid.toString(),
                   [ type: row.type, recursion: (row.recursion as Integer),
                     kmers: row.kmers.replaceAll(';', ','), coverage: (row.coverage as Integer),
                     word_size: (row.word_size?.trim() ? (row.word_size.trim() as Integer) : null) ] )
        }
        .groupTuple()

    // genetic codes for the mito annotation branch (empty if no gcode sheet -> mito annotation skipped)
    ch_gcode = params.gcode_samplesheet
        ? Channel.fromPath( params.gcode_samplesheet, checkIfExists: true )
              .splitCsv( header: true )
              .map { row -> tuple( row.taxid.toString(), (row.genetic_code as Integer) ) }
        : Channel.empty()

    SHORTREAD_ORGANELLE( ch_reads, ch_organelle )
    ORGANELLE_ANNOTATION( SHORTREAD_ORGANELLE.out.assembly, SHORTREAD_ORGANELLE.out.stats, ch_gcode )

    // ---- console summaries ----
    SHORTREAD_ORGANELLE.out.stats
        .map { meta, org, tsv -> "[assembly]      ${meta.sample}\t${org}\t" + tsv.text.trim().readLines().last() }
        .view()
    ORGANELLE_ANNOTATION.out.mito_annotation
        .map { meta, org, f -> "[mito-annot] ${meta.sample}\t${org}\t${f.name}" }.view()
    ORGANELLE_ANNOTATION.out.notes
        .map { meta, org, f -> "[note]       ${meta.sample}\t${org}\t${f.name}" }.view()
}
