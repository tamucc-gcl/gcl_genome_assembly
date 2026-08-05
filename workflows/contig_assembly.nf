/*
========================================================================================
    CONTIG_ASSEMBLY — assembler-selector subworkflow
========================================================================================
    Repo location: workflows/contig_assembly.nf

    Routes each sample to its contig assembler by meta.assembler and re-converges on a
    single per-sample assemblies channel that the downstream haplotype fork consumes
    unchanged:

      hifiasm : HiFi (+ Hi-C) -> HIFIASM   (byte-identical to the previous direct
                                            HIFIASM(ch_fastq_all) call)
      spades  : PE short reads -> SPADES    (single collapsed assembly; meta.n_hap == 1)

    The input tuple carries every contig-capable read set so each branch picks what it needs:
        tuple(meta, hifi_fastq, hic_r1, hic_r2, sr_r1, sr_r2)
    Hi-C is passed through for hifiasm phasing; sr_* for spades. Unused slots may be null
    (e.g. a hifiasm sample has null sr_r1/sr_r2).

    EMITS `assemblies` = tuple(meta, fastas), where `fastas` is the assembler's FASTA output:
      - diploid hifiasm -> a 2-element list [hap1, hap2]
      - haploid hifiasm -> a single FASTA (primary)
      - spades          -> a single FASTA (contigs)
    The fork in main.nf coerces scalar -> [x] and zips against forkHaplotypeMeta(meta),
    so all three shapes flow through untouched (spades / haploid -> one 'primary' unit).

    Branch order note: `.branch{}` sends each item to the FIRST matching selector, so the
    explicit `spades` test comes first and `hifiasm` is the catch-all (also covers an
    unset meta.assembler, defaulting to hifiasm).
========================================================================================
*/

include { HIFIASM } from '../modules/hifiasm.nf'
include { SPADES }  from '../modules/spades.nf'
include { SUBSAMPLE_SHORTREAD } from '../modules/subsample_shortread.nf'
include { COUNT_SHORTREAD } from '../modules/count_shortread.nf'

workflow CONTIG_ASSEMBLY {

    take:
    ch_reads   // tuple(meta, hifi_fastq, hic_r1, hic_r2, sr_r1, sr_r2)
    ch_traits   // tuple(sample, {telomere_motif, ploidy, haploid_genome_size})
    ch_gsize    // tuple(sample, genome_size_txt)

    main:
    ch_reads
        .branch { meta, hifi, hic1, hic2, sr1, sr2 ->
            spades:  meta.assembler == 'spades'
            hifiasm: true
        }
        .set { ch_by_assembler }

    // --- hifiasm branch: HiFi (+ optional Hi-C). Null Hi-C (HiFi-only rows) -> empty list
    //     so path staging accepts it; HIFIASM gates Hi-C phasing on meta.hic. ---
    HIFIASM(
        ch_by_assembler.hifiasm
            .map { meta, hifi, hic1, hic2, sr1, sr2 -> tuple(meta.sample, meta, hifi, hic1 ?: [], hic2 ?: []) }
            .join( ch_traits )
            .join( ch_gsize )
            .map { sample, meta, hifi, hic1, hic2, traits, gsize -> tuple(meta, hifi, hic1, hic2, traits, gsize) }
    )

    // --- spades branch: PE short reads. ---
    // genome size into meta (fallback/floor); gs file still handed to SUBSAMPLE
    ch_spades_reads = ch_by_assembler.spades
        .map { meta, hifi, h1, h2, sr1, sr2 -> tuple(meta.sample, meta, sr1, sr2) }
        .join( ch_gsize )                                      // (sample, gs_file)
        .map { s, meta, sr1, sr2, gs ->
            def g = (gs.text.trim() ==~ /\d+/) ? (gs.text.trim() as long) : 0L
            tuple(meta + [genome_size: g], sr1, sr2, gs)
        }

    SUBSAMPLE_SHORTREAD( ch_spades_reads )                     // .out.reads: (meta, r1, r2)
    COUNT_SHORTREAD( SUBSAMPLE_SHORTREAD.out.reads )           // .out.counted: (meta, r1, r2, bases.txt)

    // add measured bases; meta now carries both sizing fields for the directive
    ch_spades_in = COUNT_SHORTREAD.out.counted
        .map { meta, r1, r2, bfile ->
            tuple(meta + [sr_bases: (bfile.text.trim() as long)], r1, r2)
        }

    SPADES( ch_spades_in )

    // strip transient sizing fields -> downstream meta is pristine (no cascade)
    ch_primary = SPADES.out.contigs.map { meta, fa ->
        tuple(meta.subMap(meta.keySet() - ['genome_size', 'sr_bases']), fa)
    }

    // Re-converge on the shape the downstream fork consumes: tuple(meta, fastas).
    // (Each sample took exactly one branch, so no sample is duplicated by the mix.)
    ch_assemblies = HIFIASM.out.assemblies.mix( ch_primary )

    emit:
    assemblies = ch_assemblies
    versions   = HIFIASM.out.versions
        .mix(SPADES.out.versions)
        .mix(SUBSAMPLE_SHORTREAD.out.versions)
        .mix(COUNT_SHORTREAD.out.versions)
}
