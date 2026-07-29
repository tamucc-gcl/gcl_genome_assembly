/*
========================================================================================
    ORGANELLE — unified organelle assembly + annotation + plotting (long & short read)
========================================================================================
    Repo location: workflows/organelle.nf

    One entry point for all organelle work, branching on read type:

      HiFi  (meta.hifi)  -> MITOHIFI (assemble + annotate mito) + MITO_CIRCULAR_MAP
      other              -> SHORTREAD_ORGANELLE (GetOrganelle: mito, +plastid for plants)
                            + ORGANELLE_ANNOTATION (MITOS2 for animal/fungal mito -> GenBank
                              + circular map; GeSeq note for plant organelles)

    Input `ch_reads` = tuple(meta, hifi_fastq, sr_r1, sr_r2). Emits are unified across both
    branches as (meta, file) two-tuples (one emission per organelle), so:
      - `assemblies` feeds FILTER_ORGANELLE (grouped per sample -> bait)
      - annotation / stats / circular_map / gene_map stay compatible with the report manifest

    Everything publishes flat under ${params.outdir}/organelle with the sample in the filename.
========================================================================================
*/

include { MITOHIFI }             from '../modules/mitohifi.nf'
include { MITO_CIRCULAR_MAP }    from '../modules/mito_circular_map.nf'
include { SHORTREAD_ORGANELLE }  from './shortread_organelle.nf'
include { ORGANELLE_ANNOTATION } from './organelle_annotation.nf'

mito_circular_script = file(params.mito_circular_script ?: "${projectDir}/py_scripts/plot_mito_circular.py", checkIfExists: true)

workflow ORGANELLE {

    take:
    ch_reads       // tuple(meta, hifi_fastq, sr_r1, sr_r2)
    ch_mito_ref    // tuple(taxid, ref_fasta, ref_gb)   — MitoHiFi reference per species
    ch_gcode       // tuple(taxid, genetic_code)
    ch_organelle   // tuple(taxid, [ [type,recursion,kmers,coverage,word_size], ... ])

    main:
    ch_reads
        .branch { meta, hifi_fastq, sr1, sr2 ->
            hifi:  meta.hifi
            other: true
        }
        .set { ch_org }

    // ── HiFi branch: MitoHiFi (mito assembly + annotation) + circular map ────────────────
    ch_org.hifi
        .map { meta, hifi_fastq, sr1, sr2 -> tuple(meta.taxid?.toString(), meta, hifi_fastq) }
        .combine( ch_mito_ref, by: 0 )
        .combine( ch_gcode,    by: 0 )
        .map { taxid, meta, hifi_fastq, ref_fa, ref_gb, gcode -> tuple(meta, hifi_fastq, ref_fa, ref_gb, gcode) }
        .set { ch_mitohifi_input }
    MITOHIFI( ch_mitohifi_input )
    MITO_CIRCULAR_MAP( MITOHIFI.out.annotation, mito_circular_script )

    // ── Non-HiFi branch: GetOrganelle assembly + annotation/plotting ─────────────────────
    ch_org.other
        .map { meta, hifi_fastq, sr1, sr2 -> tuple(meta, sr1, sr2) }
        .set { ch_sr_reads }
    SHORTREAD_ORGANELLE( ch_sr_reads, ch_organelle )
    ORGANELLE_ANNOTATION( SHORTREAD_ORGANELLE.out.assembly, SHORTREAD_ORGANELLE.out.stats, ch_gcode )

    // ── Unify (2-tuples; one per organelle) ──────────────────────────────────────────────
    ch_assemblies = MITOHIFI.out.mitogenome
        .mix( SHORTREAD_ORGANELLE.out.assembly.map { meta, org, fa -> tuple(meta, fa) } )

    ch_annotation = MITOHIFI.out.annotation
        .mix( ORGANELLE_ANNOTATION.out.annotation.map { meta, org, gb -> tuple(meta, gb) } )

    ch_stats = MITOHIFI.out.stats
        .mix( ORGANELLE_ANNOTATION.out.mito_stats.map { meta, org, tsv -> tuple(meta, tsv) } )

    ch_circular = MITO_CIRCULAR_MAP.out.circular_map
        .mix( ORGANELLE_ANNOTATION.out.circular_map )

    ch_gene_map = MITOHIFI.out.gene_map
        .mix( ORGANELLE_ANNOTATION.out.gene_map.map { meta, org, png -> tuple(meta, png) } )

    emit:
    assemblies   = ch_assemblies                       // tuple(meta, fasta)  -> FILTER_ORGANELLE bait
    annotation   = ch_annotation                       // tuple(meta, gb)
    stats        = ch_stats                            // tuple(meta, tsv)    MitoHiFi-format
    circular_map = ch_circular                         // tuple(meta, png)
    gene_map     = ch_gene_map                         // tuple(meta, png)
    notes        = ORGANELLE_ANNOTATION.out.notes      // tuple(meta, org, txt)  plant organelles
    versions     = MITOHIFI.out.versions.mix( SHORTREAD_ORGANELLE.out.versions, ORGANELLE_ANNOTATION.out.versions )
}
