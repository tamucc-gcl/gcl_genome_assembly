/*
========================================================================================
    SHORTREAD_ORGANELLE — GetOrganelle assembly of organelles from short reads
========================================================================================
    Repo location: workflows/shortread_organelle.nf

    Self-contained unit: fetches the GetOrganelle DB once, fans each short-read sample out
    over its taxonomy-resolved organelle targets (mito always; + plastid for plants), and
    assembles each via GETORGANELLE (coverage-capped attempt -> optional graph prune ->
    select circle-else-linear -> status).

    Reused VERBATIM by (a) the standalone test harness (test/test_shortread_organelle.nf)
    and (b) the main pipeline's ORGANELLE_ASSEMBLY `other` branch — so behaviour is
    identical in test and production.

    take:
      ch_reads      tuple(meta, sr1, sr2)                          short-read samples
      ch_organelle  tuple(taxid, [ [type,recursion,kmers,coverage], ... ])
    emit:
      assembly  tuple(meta, organelle_type, fasta)  canonical per-organelle FASTA (absent if failed)
      graph     tuple(meta, organelle_type, gfa)    labelled assembly graph (Bandage/QC)
      stats     tuple(meta, organelle_type, tsv)    status/length/params (always, incl. failed)
      versions
========================================================================================
*/

include { DOWNLOAD_GETORGANELLE_DB } from '../modules/download_getorganelle_db.nf'
include { GETORGANELLE }            from '../modules/getorganelle.nf'

workflow SHORTREAD_ORGANELLE {
    take:
    ch_reads
    ch_organelle

    main:
    // One DB fetch, only if at least one short-read sample exists (sentinel-guarded).
    ch_db = DOWNLOAD_GETORGANELLE_DB(
                ch_reads.map { meta, sr1, sr2 -> params.getorganelle_downloads }.unique(),
                params.getorganelle_force_download
            ).config_dir.first()

    // Fan out: one GETORGANELLE task per (sample, organelle target).
    ch_reads
        .map { meta, sr1, sr2 -> tuple(meta.taxid?.toString(), meta, sr1, sr2) }
        .combine( ch_organelle, by: 0 )
        .flatMap { taxid, meta, sr1, sr2, specs -> specs.collect { s -> tuple(meta, sr1, sr2, s) } }
        .set { ch_jobs }

    GETORGANELLE( ch_jobs, ch_db )

    emit:
    assembly = GETORGANELLE.out.assembly
    graph    = GETORGANELLE.out.graph
    stats    = GETORGANELLE.out.stats
    versions = GETORGANELLE.out.versions
}
