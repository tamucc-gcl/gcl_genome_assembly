#!/usr/bin/env nextflow
/*
========================================================================================
    STANDALONE TEST HARNESS — short-read organelle assembly (SHORTREAD_ORGANELLE)
========================================================================================
    Repo location: test/test_shortread_organelle.nf

    Exercises SHORTREAD_ORGANELLE in isolation by handing it the SAME TWO channels its
    parent (ORGANELLE_ASSEMBLY, driven by main.nf) builds — nothing is re-derived here.
    This validates the real input contract:

      ch_reads      tuple(meta, sr1, sr2)                                    <- ch_reads_all (other branch)
      ch_organelle  tuple(taxid, [ [type,recursion,kmers,coverage], ... ])   <- ch_organelle_by_taxid

    In the full pipeline ch_organelle_by_taxid is computed from resolved taxonomy; here you
    PROVIDE the specs explicitly (one CSV row per organelle target per taxid), so you control
    exactly what the module receives — e.g. a single embplant_mt row for a mito-only run.
    -w / --getorganelle_from_assembly / --getorganelle_min_depth stay global params, exactly
    as in production.

    Usage:
      nextflow run test/test_shortread_organelle.nf \
          -c nextflow.config -profile slurm \
          --reads_samplesheet test/organelle_reads.csv \
          --organelle_specs   test/organelle_specs.csv \
          --outdir results_orgtest_mt

    reads sheet : sample,taxid,reads_1,reads_2
    specs sheet : taxid,type,recursion,kmers,coverage      (quote the kmers field)
========================================================================================
*/

nextflow.enable.dsl = 2

include { SHORTREAD_ORGANELLE } from '../workflows/shortread_organelle.nf'

params.reads_samplesheet = null
params.organelle_specs   = null

workflow {
    if( !params.reads_samplesheet || !params.organelle_specs )
        error "Provide --reads_samplesheet (sample,taxid,reads_1,reads_2) and --organelle_specs (taxid,type,recursion,kmers,coverage)"

    // ch_reads — exactly as ORGANELLE_ASSEMBLY hands it down: tuple(meta, sr1, sr2).
    // GETORGANELLE uses meta.sample (tag/filenames) and meta.taxid (fan-out key); a richer
    // meta from the real sample sheet would be functionally identical for this unit.
    ch_reads = Channel.fromPath( params.reads_samplesheet, checkIfExists: true )
        .splitCsv( header: true )
        .map { row ->
            tuple( [ sample: row.sample, taxid: row.taxid ],
                   file(row.reads_1, checkIfExists: true),
                   file(row.reads_2, checkIfExists: true) )
        }

    // ch_organelle — exactly as ch_organelle_by_taxid: tuple(taxid, [ spec-map, ... ]).
    // One CSV row per (taxid, organelle target); grouped into the per-taxid spec list.
    ch_organelle = Channel.fromPath( params.organelle_specs, checkIfExists: true )
        .splitCsv( header: true )
        .map { row ->
            tuple( row.taxid.toString(),
                   [ type      : row.type,
                     recursion : (row.recursion as Integer),
                     kmers     : row.kmers,
                     coverage  : (row.coverage as Integer) ] )
        }
        .groupTuple()

    SHORTREAD_ORGANELLE( ch_reads, ch_organelle )

    // console summary of the status contract (also published under --outdir)
    SHORTREAD_ORGANELLE.out.stats
        .map { meta, org_type, tsv -> "[status] ${meta.sample}\t${org_type}\t" + tsv.text.trim().readLines().last() }
        .view()
}
