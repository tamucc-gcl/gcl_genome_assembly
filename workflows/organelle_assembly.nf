/*
========================================================================================
    ORGANELLE_ASSEMBLY — organelle assembly + annotation (selector subworkflow)
========================================================================================
    Repo location: workflows/organelle_assembly.nf

    Abstracts organelle handling so the contig path routes through one place regardless of
    read type. Branches on meta.hifi:

      HiFi  -> MitoHiFi (mito only) + circular map
      other -> SHORTREAD_ORGANELLE (GetOrganelle): mito for animals/fungi, plastid + mito for
               plants, via the taxonomy-resolved organelle-spec side-channel.

    Input `ch_reads` = tuple(meta, hifi_fastq, sr_r1, sr_r2) — carries both read sets so each
    branch picks what it needs (HiFi fastq for MitoHiFi; sr_* for GetOrganelle).

    EMITS the MitoHiFi surface main.nf already consumes, PLUS the short-read organelle outputs:
      mitogenome / annotation / stats / circular_map          (HiFi/MitoHiFi — unchanged)
      sr_organelle / sr_organelle_graph / sr_organelle_stats  (GetOrganelle)
========================================================================================
*/

include { MITOHIFI }            from '../modules/mitohifi.nf'
include { MITO_CIRCULAR_MAP }   from '../modules/mito_circular_map.nf'
include { SHORTREAD_ORGANELLE } from './shortread_organelle.nf'
mito_circular_script = file("${projectDir}/py_scripts/plot_mito_circular.py", checkIfExists: true)

workflow ORGANELLE_ASSEMBLY {

    take:
    ch_reads       // tuple(meta, hifi_fastq, sr_r1, sr_r2)
    ch_mito_ref    // tuple(taxid, ref_fasta, ref_gb) — per resolved species
    ch_gcode       // tuple(taxid, genetic_code)
    ch_organelle   // tuple(taxid, [ [type,recursion,kmers,coverage], ... ])

    main:
    ch_reads
        .branch { meta, hifi_fastq, sr1, sr2 ->
            hifi:  meta.hifi
            other: true
        }
        .set { ch_org }

    // --- HiFi branch: MitoHiFi (mito only), each sample paired with its own species' reference.
    ch_org.hifi
        .map { meta, hifi_fastq, sr1, sr2 -> tuple(meta.taxid?.toString(), meta, hifi_fastq) }
        .combine( ch_mito_ref, by: 0 )
        .combine( ch_gcode,    by: 0 )
        .map { taxid, meta, hifi_fastq, ref_fa, ref_gb, gcode -> tuple(meta, hifi_fastq, ref_fa, ref_gb, gcode) }
        .set { ch_mitohifi_input }
    MITOHIFI(ch_mitohifi_input)

    MITO_CIRCULAR_MAP(MITOHIFI.out.annotation, mito_circular_script)

    // --- Non-HiFi branch: GetOrganelle on short reads (mito always; + plastid for plants). ---
    ch_org.other
        .map { meta, hifi_fastq, sr1, sr2 -> tuple(meta, sr1, sr2) }
        .set { ch_sr_reads }
    SHORTREAD_ORGANELLE( ch_sr_reads, ch_organelle )

    emit:
    mitogenome   = MITOHIFI.out.mitogenome
    annotation   = MITOHIFI.out.annotation
    stats        = MITOHIFI.out.stats
    circular_map = MITO_CIRCULAR_MAP.out.circular_map

    sr_organelle       = SHORTREAD_ORGANELLE.out.assembly
    sr_organelle_graph = SHORTREAD_ORGANELLE.out.graph
    sr_organelle_stats = SHORTREAD_ORGANELLE.out.stats

    versions = MITOHIFI.out.versions.mix( SHORTREAD_ORGANELLE.out.versions )
}