/*
========================================================================================
    QC_PHASE SUBWORKFLOW
========================================================================================
    Repo location: workflows/qc_phase.nf

    Per-stage assembly QC across every checkpoint, plus the cross-stage compilation and the
    interactive HTML viewer.

    LEVEL-2 INTERFACE
      Callers pass ONE labeled channel `staged_assemblies` of (meta, stage, fasta) — each
      element is one haplotype's assembly at one stage. The stage key routes each element to
      its ASSEMBLY_QC alias via .branch. The aliases and the ASSEMBLY_QC subworkflow are
      UNCHANGED, so per-stage task caching (BUSCO/QUAST/MERQURY/MAPPING) is preserved.

      Stage guards (run_all_qc / run_<x>) live at the CALLER, in how `staged_assemblies` is
      built: a stage is QC'd iff the caller mixed it in. 'final' is always present. A stage
      that isn't mixed in arrives as an empty branch -> its ASSEMBLY_QC alias runs zero
      tasks (== not called), so behaviour matches the old per-call `if` guards exactly.

      Hi-C bam/pairs metrics are produced upstream (scaffolding phase) and passed in for the
      compile step.

    Stage keys (must match how the caller tags `staged_assemblies`):
      initial, organelle_filtered, purged, redundans, contig_corrected, contig_decontam,
      scaffold, scaffold_corrected, scaffold_decontam, scaffold_round2, gap_filled,
      teloclip, final
========================================================================================
*/

nextflow.enable.dsl = 2

include { ASSEMBLY_QC as ASSEMBLY_QC_INITIAL            } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_ORGANELLE_FILTERED } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_PURGED             } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_REDUNDANS          } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_CONTIG_CORRECTED   } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_CONTIG_DECONTAM    } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD           } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_CORRECTED } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_DECONTAM  } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_SCAFFOLD_ROUND2    } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_GAP_FILLED         } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_TELOCLIP           } from './assembly_qc.nf'
include { ASSEMBLY_QC as ASSEMBLY_QC_FINAL              } from './assembly_qc.nf'
include { COMPILE_FINAL_QC } from '../modules/compile_final_qc.nf'
include { ASSEMBLY_REPORT  } from '../modules/assemblyReport.nf'

workflow QC_PHASE {
    take:
    staged_assemblies       // (meta, stage, fasta)         — per-haplotype, per-stage; guards applied by caller
    hifi_reads              // (meta, hifi_fastq)            — per-sample
    meryl_db                // (meta, meryl_db)              — per-sample
    busco_db                // value: MAP taxid -> busco_lineage
    all_bam_metrics         // (meta, checkpoint, tsv)       — Hi-C BAM metrics from scaffolding phase
    all_pairs_metrics       // (meta, checkpoint, tsv)       — Hi-C pairs metrics from scaffolding phase
    compile_qc_script       // value: compile_qc.R
    assembly_report_script  // value: generate_assembly_report.py

    main:

    // Route each stage to its ASSEMBLY_QC alias (stage key is unique per checkpoint).
    staged_assemblies
        .branch {
            initial:            it[1] == 'initial'
            organelle_filtered: it[1] == 'organelle_filtered'
            purged:             it[1] == 'purged'
            redundans:          it[1] == 'redundans'
            contig_corrected:   it[1] == 'contig_corrected'
            contig_decontam:    it[1] == 'contig_decontam'
            scaffold:           it[1] == 'scaffold'
            scaffold_corrected: it[1] == 'scaffold_corrected'
            scaffold_decontam:  it[1] == 'scaffold_decontam'
            scaffold_round2:    it[1] == 'scaffold_round2'
            gap_filled:         it[1] == 'gap_filled'
            teloclip:           it[1] == 'teloclip'
            finalstage:         it[1] == 'final'
        }
        .set { st }

    ch_all_assembly_summaries = Channel.empty()

    ASSEMBLY_QC_INITIAL(            st.initial.map            { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'contig')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_INITIAL.out.assembly_summary)

    ASSEMBLY_QC_ORGANELLE_FILTERED(st.organelle_filtered.map { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'contig_organelle_filtered')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_ORGANELLE_FILTERED.out.assembly_summary)

    ASSEMBLY_QC_PURGED(            st.purged.map             { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'contig_purged')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_PURGED.out.assembly_summary)

    ASSEMBLY_QC_REDUNDANS(         st.redundans.map          { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'contig_purged')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_REDUNDANS.out.assembly_summary)

    ASSEMBLY_QC_CONTIG_CORRECTED(  st.contig_corrected.map   { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'contig_corrected')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_CONTIG_CORRECTED.out.assembly_summary)

    ASSEMBLY_QC_CONTIG_DECONTAM(   st.contig_decontam.map    { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'contig_decontam')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_CONTIG_DECONTAM.out.assembly_summary)

    ASSEMBLY_QC_SCAFFOLD(          st.scaffold.map           { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'scaffold')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_SCAFFOLD.out.assembly_summary)

    ASSEMBLY_QC_SCAFFOLD_CORRECTED(st.scaffold_corrected.map { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'scaffold_corrected')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_SCAFFOLD_CORRECTED.out.assembly_summary)

    ASSEMBLY_QC_SCAFFOLD_DECONTAM( st.scaffold_decontam.map  { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'scaffold_decontam')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_SCAFFOLD_DECONTAM.out.assembly_summary)

    ASSEMBLY_QC_SCAFFOLD_ROUND2(   st.scaffold_round2.map    { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'scaffold_round2')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_SCAFFOLD_ROUND2.out.assembly_summary)

    ASSEMBLY_QC_GAP_FILLED(        st.gap_filled.map         { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'gap_filled')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_GAP_FILLED.out.assembly_summary)

    ASSEMBLY_QC_TELOCLIP(          st.teloclip.map           { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'teloclip_extended')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_TELOCLIP.out.assembly_summary)

    ASSEMBLY_QC_FINAL(             st.finalstage.map         { m, s, f -> tuple(m, f) }, hifi_reads, meryl_db, busco_db, 'final')
    ch_all_assembly_summaries = ch_all_assembly_summaries.mix(ASSEMBLY_QC_FINAL.out.assembly_summary)

    ch_final_busco = ASSEMBLY_QC_FINAL.out.busco_results
    ch_final_busco_table = ASSEMBLY_QC_FINAL.out.busco_full_table

    // Compile cross-stage report: summaries built here + Hi-C metrics from the scaffolding phase.
    COMPILE_FINAL_QC(
        ch_all_assembly_summaries.map { sample_id, qc_label, tsv -> tsv }.collect().ifEmpty([]),
        all_bam_metrics.map           { meta, checkpoint, tsv -> tsv }.collect().ifEmpty([]),
        all_pairs_metrics.map         { meta, checkpoint, tsv -> tsv }.collect().ifEmpty([]),
        compile_qc_script
    )

    // Interactive HTML assembly viewer.
    ASSEMBLY_REPORT(COMPILE_FINAL_QC.out.metrics, assembly_report_script)

    emit:
    metrics     = COMPILE_FINAL_QC.out.metrics       // assembly_qc_metrics.csv
    plots       = COMPILE_FINAL_QC.out.plots         // *.png
    report_html = ASSEMBLY_REPORT.out.report_html    // assembly_qc_report.html
    final_busco = ch_final_busco                     // per-haplotype (meta, results) — for SNAIL
    final_busco_table = ch_final_busco_table    // per-haplotype (meta, full_table.tsv) — for SNAIL
    versions    = ASSEMBLY_QC_FINAL.out.versions
}
