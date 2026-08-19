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

mito_circular_script = file("${projectDir}/py_scripts/plot_mito_circular.py", checkIfExists: true)

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

    // ---- Per-sample bait bundles for FILTER_ORGANELLE ---------------------------------
    // Exactly one emission per input sample (empty list = produced no organelle), with the
    // gather confined to a SINGLE read-type branch. main.nf used to gather globally over
    // ch_assemblies, so one slow or newly added sample in either branch stalled contig
    // filtering -- and every step downstream of it -- for every other sample.
    //
    // Both arms end in groupTuple() so the bait value stays the same ArrayBag the global
    // gather produced. A bare [fa] list would be a different collection type and is not
    // guaranteed to hash identically, which would re-run FILTER_ORGANELLE and cascade from
    // there through every downstream task.
    //
    //   HiFi       MITOHIFI.out.mitogenome is NOT optional and errorStrategy is
    //              retry/terminate (never 'ignore'), so a HiFi sample yields exactly one
    //              mitogenome or the run fails. groupKey(sample, 1) therefore lets
    //              groupTuple emit on that single arrival with no dependence on any channel
    //              closing. If MITOHIFI ever gains 'ignore' or an optional mitogenome this
    //              arm needs the remainder join the short-read arm uses, or a sample will
    //              go missing from the combine in main.nf.
    //   short-read GETORGANELLE.out.assembly IS optional and a sample can have >1 organelle
    //              target, so the count is not knowable up front and this arm gathers --
    //              but only over ch_org.other. The remainder join re-attaches an empty list
    //              for samples that produced none. Sorting is applied only above size 1,
    //              where groupTuple's arrival order is genuinely unstable; at size 1 it
    //              would convert the ArrayBag to a List for no benefit.
    ch_hifi_baits = MITOHIFI.out.mitogenome
        .map { meta, fa -> tuple(groupKey(meta.sample, 1), meta.sample, fa) }
        .groupTuple()
        .map { key, samples, fas -> tuple(samples[0], fas) }

    ch_sr_baits = ch_org.other
        .map { meta, hifi_fastq, sr1, sr2 -> tuple(meta.sample, meta.sample) }
        .join( SHORTREAD_ORGANELLE.out.assembly
                   .map { meta, org, fa -> tuple(meta.sample, fa) }
                   .groupTuple(),
               remainder: true )
        .map { sample, s, fas ->
            tuple(sample, (fas == null) ? [] : (fas.size() > 1 ? fas.sort { it.name } : fas)) }

    ch_baits = ch_hifi_baits.mix( ch_sr_baits )

    emit:
    assemblies   = ch_assemblies                       // tuple(meta, fasta)  -> FILTER_ORGANELLE bait
    baits        = ch_baits                            // tuple(sample, [fasta,...]) per sample
    annotation   = ch_annotation                       // tuple(meta, gb)
    stats        = ch_stats                            // tuple(meta, tsv)    MitoHiFi-format
    circular_map = ch_circular                         // tuple(meta, png)
    gene_map     = ch_gene_map                         // tuple(meta, png)
    notes        = ORGANELLE_ANNOTATION.out.notes      // tuple(meta, org, txt)  plant organelles
    versions     = MITOHIFI.out.versions.mix( SHORTREAD_ORGANELLE.out.versions, ORGANELLE_ANNOTATION.out.versions )
}
