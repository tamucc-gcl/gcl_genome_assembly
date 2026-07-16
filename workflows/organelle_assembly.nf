/*
========================================================================================
    ORGANELLE_ASSEMBLY — organelle assembly + annotation (selector subworkflow)
========================================================================================
    Repo location: workflows/organelle_assembly.nf

    Abstracts organelle handling so the contig path routes through one place regardless of
    read type. Branches on meta.hifi:

      HiFi  -> MitoHiFi (mito only) + circular map        [behaviour identical to the old
                                                            STEP 3b direct calls]
      other -> STUB (Phase 4b: GetOrganelle on short reads — mito for animals, mito +
               chloroplast for plants via the species/taxid kingdom lookup)

    4a-ii scope: the non-HiFi branch is a NO-OP placeholder (not consumed), so short-read
    samples currently carry no organelle and are dropped at FILTER_MITO_CONTIGS downstream
    (their contig path is gated/branched in 4a-iii). HiFi behaviour is unchanged. 4b fills
    the `other` branch and (likely) moves the organelle-contig filtering in here too.

    Input `ch_reads` = tuple(meta, hifi_fastq, sr_r1, sr_r2) — carries both read sets so each
    branch can pick what it needs (HiFi fastq for MitoHiFi; sr_* for GetOrganelle in 4b).
    Reference channels come from the (currently global) FIND_MITO_REFERENCE.

    EMITS the same surface main.nf consumed from MITOHIFI/MITO_CIRCULAR_MAP:
      mitogenome   -> FILTER_MITO_CONTIGS
      annotation   -> manifest (mito .gb)
      stats        -> manifest + report
      circular_map -> manifest (mito gene-map image)
    (empty for a run with no HiFi samples — downstream consumers simply receive nothing.)
========================================================================================
*/

include { MITOHIFI }          from '../modules/mitohifi.nf'
include { MITO_CIRCULAR_MAP } from '../modules/mito_circular_map.nf'
mito_circular_script = file("${projectDir}/py_scripts/plot_mito_circular.py",           checkIfExists: true)

workflow ORGANELLE_ASSEMBLY {

    take:
    ch_reads       // tuple(meta, hifi_fastq, sr_r1, sr_r2)
    ch_mito_ref    // tuple(taxid, ref_fasta, ref_gb) — per resolved species

    main:
    ch_reads
        .branch { meta, hifi_fastq, sr1, sr2 ->
            hifi:  meta.hifi
            other: true
        }
        .set { ch_org }

    // --- HiFi branch: MitoHiFi (mito only), each sample paired with its own species' reference.
    // combine(by:0) is one-to-many (both haplotypes of a diploid sample share the one reference).
    // A HiFi sample whose taxid resolved no reference is dropped here (carries no mito) — same as
    // the existing behaviour for organelle-less samples at FILTER_MITO_CONTIGS.
    ch_org.hifi
        .map { meta, hifi_fastq, sr1, sr2 -> tuple(meta.taxid?.toString(), meta, hifi_fastq) }
        .combine( ch_mito_ref, by: 0 )
        .map { taxid, meta, hifi_fastq, ref_fa, ref_gb -> tuple(meta, hifi_fastq, ref_fa, ref_gb) }
        .set { ch_mitohifi_input }
    MITOHIFI(ch_mitohifi_input)

    MITO_CIRCULAR_MAP(MITOHIFI.out.annotation, mito_circular_script)

    // --- Non-HiFi branch: STUB (filled in Phase 4b with GetOrganelle). ---
    // Intentionally not consumed yet; short-read samples carry no organelle in 4a-ii.
    // ch_org.other

    emit:
    mitogenome   = MITOHIFI.out.mitogenome
    annotation   = MITOHIFI.out.annotation
    stats        = MITOHIFI.out.stats
    circular_map = MITO_CIRCULAR_MAP.out.circular_map
}
