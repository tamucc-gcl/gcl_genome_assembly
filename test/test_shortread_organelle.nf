#!/usr/bin/env nextflow
/*
========================================================================================
    STANDALONE TEST HARNESS — short-read organelle assembly
========================================================================================
    Repo location: test/test_shortread_organelle.nf

    Runs SHORTREAD_ORGANELLE in isolation (no BAM->fastq, no hifiasm, no QC), so the
    GetOrganelle path can be developed and tested on its own, then dropped into the main
    pipeline unchanged — main.nf's ORGANELLE_ASSEMBLY calls the SAME subworkflow.

    Usage:
      nextflow run test/test_shortread_organelle.nf \\
          -c nextflow.config \\
          --organelle_samplesheet test/organelle_samplesheet.csv \\
          --outdir results_orgtest
      # add -stub for a wiring smoke-test (no real assembly)

    Samplesheet (header required):
      sample,taxid,kingdom,reads_1,reads_2
        kingdom : plant | animal | fungi | other   (plant => pt+mt targets; else mt)
        taxid   : grouping key only in this harness (taxonomy is NOT resolved here); set it
                  to the real NCBI taxid so it matches how the full pipeline keys things
        reads_* : trimmed Illumina PE fastq(.gz)
========================================================================================
*/

nextflow.enable.dsl = 2

include { SHORTREAD_ORGANELLE } from '../workflows/shortread_organelle.nf'
include { organelleTypesFor; getorganelleRecursionFor; getorganelleKmersFor; getorganelleCoverageFor } from '../functions/taxonomy.nf'

params.organelle_samplesheet = null

workflow {
    if( !params.organelle_samplesheet )
        error "Provide --organelle_samplesheet <csv> with columns: sample,taxid,kingdom,reads_1,reads_2"

    // short-read samples: tuple(meta, r1, r2)
    ch_reads = Channel.fromPath( params.organelle_samplesheet, checkIfExists: true )
        .splitCsv( header: true )
        .map { row ->
            tuple( [ sample: row.sample, taxid: row.taxid ],
                   file(row.reads_1, checkIfExists: true),
                   file(row.reads_2, checkIfExists: true) )
        }

    // per-taxid organelle specs, built from the SAME helpers the pipeline uses (kingdom-driven
    // here). Read the sheet a second time so we never double-consume a single channel.
    ch_organelle = Channel.fromPath( params.organelle_samplesheet, checkIfExists: true )
        .splitCsv( header: true )
        .map { row ->
            def types = params.getorganelle_organelle_types
                            ? params.getorganelle_organelle_types.tokenize(',')*.trim()
                            : organelleTypesFor( [ kingdom: row.kingdom ] )
            def specs = types.collect { t ->
                [ type      : t,
                  recursion : (params.getorganelle_recursion ?: getorganelleRecursionFor(t)),
                  kmers     : (params.getorganelle_kmers     ?: getorganelleKmersFor(t)),
                  coverage  : (params.getorganelle_coverage  ?: getorganelleCoverageFor(t)) ]
            }
            tuple( row.taxid?.toString(), specs )
        }
        .unique()

    SHORTREAD_ORGANELLE( ch_reads, ch_organelle )

    // console summary of the status contract (also published to --outdir)
    SHORTREAD_ORGANELLE.out.stats
        .map { meta, org_type, tsv -> "[status] ${meta.sample}\t${org_type}\t" + tsv.text.trim().readLines().last() }
        .view()
}
