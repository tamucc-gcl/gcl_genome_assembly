#!/usr/bin/env nextflow
// Standalone test for SHORTREAD_ORGANELLE — provides ch_reads + ch_organelle explicitly.
//   reads sheet: sample,taxid,reads_1,reads_2
//   specs sheet: taxid,type,recursion,kmers,coverage,word_size   (k-list uses ';'; word_size blank => auto)
nextflow.enable.dsl = 2

include { SHORTREAD_ORGANELLE } from '../workflows/shortread_organelle.nf'

params.reads_samplesheet = null
params.organelle_specs   = null

workflow {
    if( !params.reads_samplesheet || !params.organelle_specs )
        error "Provide --reads_samplesheet (sample,taxid,reads_1,reads_2) and --organelle_specs (taxid,type,recursion,kmers,coverage,word_size)"

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

    SHORTREAD_ORGANELLE( ch_reads, ch_organelle )

    SHORTREAD_ORGANELLE.out.stats
        .map { meta, org_type, tsv -> "[status] ${meta.sample}\t${org_type}\t" + tsv.text.trim().readLines().last() }
        .view()
}
