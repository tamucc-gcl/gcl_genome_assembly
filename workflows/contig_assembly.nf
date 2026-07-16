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

workflow CONTIG_ASSEMBLY {

    take:
    ch_reads   // tuple(meta, hifi_fastq, hic_r1, hic_r2, sr_r1, sr_r2)

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
        ch_by_assembler.hifiasm.map { meta, hifi, hic1, hic2, sr1, sr2 -> tuple(meta, hifi, hic1 ?: [], hic2 ?: []) }
    )

    // --- spades branch: PE short reads. ---
    SPADES(
        ch_by_assembler.spades.map { meta, hifi, hic1, hic2, sr1, sr2 -> tuple(meta, sr1, sr2) }
    )

    // Re-converge on the shape the downstream fork consumes: tuple(meta, fastas).
    // (Each sample took exactly one branch, so no sample is duplicated by the mix.)
    ch_assemblies = HIFIASM.out.assemblies.mix( SPADES.out.contigs )

    emit:
    assemblies = ch_assemblies
}
