/*
========================================================================================
    REPORTING  —  report manifest assembly + summary report (terminal pipeline tail)
========================================================================================
    Repo location:  workflows/reporting.nf

    Extracted verbatim from the tail of the main `workflow {}` (the "SUMMARY REPORT —
    Build manifest" section through the SUMMARY_REPORT call) so the main workflow body
    stays under Groovy's 65,535-char compiled-constant limit (Nextflow captures each
    workflow body's source text as a String constant).

    Purely terminal: it consumes already-produced upstream channels, builds the report
    manifest + provenance TSVs, and calls SUMMARY_REPORT. It emits nothing.

    Every `take:` below is a channel that was a bare `ch_*` variable or a `PROCESS.out.*`
    reference in main.nf; the process that produces it stays invoked in main and its
    output is passed in here. The `take:` order is the call order — keep them in sync.
========================================================================================
*/

include { COLLECT_SOFTWARE_VERSIONS } from '../modules/collect_software_versions.nf'
include { SUMMARY_REPORT           } from '../modules/summary_report'

// ---------------------------------------------------------------------------
// Parse heterozygosity (%) and repeat content (%) from a GenomeScope2
// summary.txt. GenomeScope prints heterozygosity as a min-max range and does
// NOT print a repeat percentage, so repeat content is derived as
//     Genome Repeat Length / Genome Haploid Length * 100
// Both take the MAX column, matching ESTIMATE_GENOME_SIZE's genome-size parse
// (awk '{print $(NF-1)}'). Returns [het, rep] as ready-to-print % strings, or
// 'NA' when a field is absent (higher-ploidy models, stub runs, failed fits).
// File scope (not the workflow body) so it does not grow the workflow's
// compiled source-String toward Groovy's 65,535-char limit.
// ---------------------------------------------------------------------------
def parseGenomescopeSummary(gs_summary) {
    def out = [het: 'NA', rep: 'NA']
    def summary = file("${gs_summary}")
    if (!summary.exists()) return out

    Long hapLen = null
    Long repLen = null
    try {
        summary.readLines().each { String line ->
            def mh = (line =~ /Heterozygous \(ab\)\s+[0-9.]+%\s+([0-9.]+)%/)
            if (mh.find()) out.het = String.format(Locale.US, '%.2f%%', (mh.group(1) as double))

            def mg = (line =~ /Genome Haploid Length\s+[0-9,]+\s*bp\s+([0-9,]+)\s*bp/)
            if (mg.find()) hapLen = (mg.group(1).replaceAll(',', '') as Long)

            def mr = (line =~ /Genome Repeat Length\s+[0-9,]+\s*bp\s+([0-9,]+)\s*bp/)
            if (mr.find()) repLen = (mr.group(1).replaceAll(',', '') as Long)
        }
    } catch (Exception e) {
        return out
    }
    if (repLen != null && hapLen != null && hapLen > 0L)
        out.rep = String.format(Locale.US, '%.1f%%', (repLen / (double) hapLen) * 100.0d)
    return out
}

workflow REPORTING {

    take:
    ch_finalized_assembly          //  FINALIZE_ASSEMBLY.out.assembly     tuple(meta, fasta)
    ch_snail                       //  SNAIL_PLOT_FINAL.out.snail          tuple(id, qc_label, svg)
    ch_contact_maps                //  CONTACT_MAP_FINAL.out.contact_maps  tuple(meta, stage, pngs) | empty
    ch_dotplot                     //  FINAL_VIZ.out.dotplot               tuple(id1, id2, png)
    ch_riparian                    //  FINAL_VIZ.out.riparian              tuple(id1, id2, png)
    ch_tidk_plot                   //  FINAL_VIZ.out.tidk_plot             tuple(id, svg)
    ch_compiled_metrics            //  COMPILE_FINAL_QC.out.metrics        path (also SUMMARY_REPORT arg 2)
    ch_compiled_plots              //  COMPILE_FINAL_QC.out.plots          path(s)
    ch_assembly_report_html        //  ASSEMBLY_REPORT.out.report_html     path
    ch_organelle_annotation        //  ORGANELLE.out.annotation            tuple(meta, gb)
    ch_organelle_stats             //  ORGANELLE.out.stats                 tuple(meta, tsv) (used twice)
    ch_organelle_circular          //  ORGANELLE.out.circular_map          tuple(meta, png)
    ch_genomescope_results         //  ESTIMATE_GENOME_SIZE.out.results    tuple(meta, summary_txt)
    ch_genome_size_est             //  ESTIMATE_GENOME_SIZE.out.size       tuple(meta, size_file)
    ch_sample_identity             //  per-sample taxonomy side-channel    tuple(sample, tax)
    ch_input                       //  parsed sample sheet                 tuple(meta, reads)
    ch_ploidy_by_sample            //  per-sample ploidy side-channel      tuple(sample, ploidy)
    ch_telomere_for_report         //  FINAL_VIZ.out.telomere_summary
    ch_pairwise_summary            //  FINAL_VIZ.out.pairwise_summary
    ch_teloclip_stats_for_report   //  COLLECT_TELOCLIP_STATS.out.stats | NO_TELOCLIP
    ch_pangenome_report            //  PANGENOME.out.report | NO_PANGENOME
    ch_name_map                    //  COLLECT_NAME_MAPS.out.map | NO_NAMEMAP
    ch_versions                    //  accumulated versions channel (already fully mixed in main)
    ch_summary_report_script       //  file: r_scripts/generate_summary_report.R

    main:

    // =========================================================================
    //  SUMMARY REPORT — Build manifest and call the process
    //
    //  Verified publishDir paths from each module:
    //    FINALIZE_ASSEMBLY  → ${params.outdir}/assembly/final
    //    GAP_FILLING        → ${params.outdir}/assembly/scaffold/gap_filling
    //    SNAIL_PLOT         → ${params.outdir}/snail_plots
    //    CONTACT_MAP        → ${params.outdir}/contact_maps
    //    PAIRWISE_ALIGNMENT → ${params.outdir}/pairwise_alignments
    //    COMPILE_FINAL_QC   → ${params.outdir}/qc/assembly
    //    ASSEMBLY_REPORT    → ${params.outdir}/reports
    // =========================================================================

    // ---- Final genome assemblies ----
    // FINALIZE_ASSEMBLY.out.assembly: tuple(haplotype_id, fasta)
    // publishDir: ${params.outdir}/assembly/final
    ch_manifest_assemblies = ch_finalized_assembly
        .map { meta, fasta -> "assembly\t${meta.id}\t.\t${fasta.name}\tassembly/final" }

    // ---- Snail plots (final) ----
    // SNAIL_PLOT_FINAL.out.snail: tuple(haplotype_id, qc_label, svg)
    // publishDir: ${params.outdir}/snail_plots
    ch_manifest_snails = ch_snail
        .map { hap_id, qc_label, svg ->
            "snail\t${hap_id}\t.\t${svg.name}\tsnail_plots"
        }

    // ---- Contact maps (conditional) ----
    // CONTACT_MAP_FINAL.out.contact_maps: tuple(haplotype_id, stage, png_files)
    // publishDir: ${params.outdir}/contact_maps
    if (params.run_final_contact_maps) {
        ch_manifest_contact_maps = ch_contact_maps
            .flatMap { meta, stage, pngs ->
                def png_list = pngs instanceof List ? pngs : [pngs]
                png_list.collect { png ->
                    "contact_map\t${meta.id}\t.\t${png.name}\tcontact_maps"
                }
            }
    } else {
        ch_manifest_contact_maps = Channel.empty()
    }

    // ---- Pairwise dotplots (conditional) ----
    // FINAL_VIZ.out.dotplot: tuple(id1, id2, png)
    // publishDir: ${params.outdir}/pairwise_alignments
    if (params.run_pairwise_alignments) {
        ch_manifest_dotplots = ch_dotplot
            .map { id1, id2, png ->
                "dotplot\t${id1}\t${id2}\t${png.name}\tpairwise_alignments"
            }
    } else {
        ch_manifest_dotplots = Channel.empty()
    }

    // ---- Pairwise riparian plots (conditional) ----
    // FINAL_VIZ.out.riparian: tuple(id1, id2, png)
    // publishDir: ${params.outdir}/pairwise_alignments
    if (params.run_pairwise_alignments) {
        ch_manifest_riparian = ch_riparian
            .map { id1, id2, png ->
                "riparian\t${id1}\t${id2}\t${png.name}\tpairwise_alignments"
            }
    } else {
        ch_manifest_riparian = Channel.empty()
    }

    // ---- tidk plot SVGs for manifest ----
    ch_manifest_tidk_plots = ch_tidk_plot
        .map { hap_id, svg ->
            "tidk_plot\t${hap_id}\t.\t${svg.name}\ttelomeres/plots"
        }


    // ---- Compiled QC CSV from COMPILE_FINAL_QC ----
    // publishDir: ${params.outdir}/qc/assembly
    ch_manifest_compiled_csv = ch_compiled_metrics
        .map { csv ->
            "compiled_qc\t.\t.\t${csv.name}\tqc/assembly"
        }

    // ---- QC trend plots (PNGs) from COMPILE_FINAL_QC ----
    // publishDir: ${params.outdir}/qc/assembly
    ch_manifest_qc_plots = ch_compiled_plots
        .flatten()
        .map { png ->
            "qc_plot\t.\t.\t${png.name}\tqc/assembly"
        }

    // ---- Interactive HTML report from ASSEMBLY_REPORT ----
    // publishDir: ${params.outdir}/reports
    ch_manifest_assembly_report = ch_assembly_report_html
        .map { html ->
            "assembly_report_html\t.\t.\t${html.name}\treports"
        }

    // ---- Mitogenome assembly ----
    // publishDir: ${params.outdir}/mitogenome/${sample_id}
    ch_manifest_mito_gb = ch_organelle_annotation
        .map { meta, gb -> "mito_genbank\t${meta.sample}\t.\t${gb.name}\torganelle" }

    ch_manifest_mito_stats = ch_organelle_stats
        .map { meta, tsv -> "mito_stats\t${meta.sample}\t.\t${tsv.name}\torganelle" }

    ch_manifest_mito_circular = ch_organelle_circular
        .map { meta, png -> "mito_gene_map\t${meta.sample}\t.\t${png.name}\torganelle" }

    // ---- GenomeScope profiles (linear_plot.png per sample) ----
    // ESTIMATE_GENOME_SIZE.out.results: tuple(meta, <sample>.summary.txt)
    // publishDir: ${params.outdir}/est_genome_size (flat, ID-prefixed) → est_genome_size/<sample>.linear_plot.png
    ch_manifest_genomescope = ch_genomescope_results
        .map { meta, gs_summary ->
            "genomescope_plot\t${meta.sample}\t.\t${meta.sample}.linear_plot.png\test_genome_size"
        }

    // ---- Combine all manifest entries into a single TSV ----
    ch_manifest_assemblies
        .mix(ch_manifest_snails)
        .mix(ch_manifest_contact_maps)
        .mix(ch_manifest_dotplots)
        .mix(ch_manifest_riparian)
        .mix(ch_manifest_compiled_csv)
        .mix(ch_manifest_qc_plots)
        .mix(ch_manifest_genomescope)
        .mix(ch_manifest_assembly_report)
        .mix(ch_manifest_mito_gb)
        .mix(ch_manifest_mito_stats)
        .mix(ch_manifest_mito_circular)
        .mix(ch_manifest_tidk_plots)
        .collectFile(
            name: 'report_manifest.tsv',
            seed: 'type\tid\tid2\tfilename\tsubdir',
            newLine: true
        )
        .set { ch_report_manifest }

    // Collect mitogenome stats for report
    ch_mito_stats_for_report = ch_organelle_stats
        .map { meta, tsv -> tsv }
        .collectFile(
            name: 'all_mito_stats.tsv',
            keepHeader: true,
            skip: 1
        )
        .ifEmpty(file('NO_MITO_STATS'))

    // ---- Per-sample taxonomy + genome-size estimate (4b-i Increment 4) ----
    ch_sample_taxonomy_tsv = ch_sample_identity
        .map { sample, tax ->
            "${sample}\t${tax.taxid}\t${tax.name}\t${tax.kingdom}\t${tax.busco_lineage}\t${tax.genetic_code}\t${tax.telomere_motif}" }
        .collectFile(name: 'sample_taxonomy.tsv', newLine: true,
                     seed: 'sample\ttaxid\tspecies\tkingdom\tbusco_lineage\tgenetic_code\ttelomere_motif',
                     sort: true)
        .ifEmpty(file('NO_TAXONOMY'))

    ch_genome_size_tsv = ch_genome_size_est
        .map { meta, size_file -> "${meta.sample}\t${size_file.text.trim()}" }
        .collectFile(name: 'genome_sizes.tsv',
                     seed: 'sample\test_genome_size_bp',
                     newLine: true)
        .ifEmpty(file('NO_GENOME_SIZE'))

    // ---- Per-sample GenomeScope heterozygosity + repeat content ----
    // Parsed straight from each sample's genomescope summary.txt (already flowing
    // as ch_genomescope_results) via parseGenomescopeSummary() defined above, so
    // nothing changes in ESTIMATE_GENOME_SIZE or main.nf. Mirrors genome_sizes.tsv.
    ch_genomescope_metrics_tsv = ch_genomescope_results
        .map { meta, gs_summary ->
            def gm = parseGenomescopeSummary(gs_summary)
            "${meta.sample}\t${gm.het}\t${gm.rep}"
        }
        .collectFile(name: 'genomescope_metrics.tsv',
                     seed: 'sample\theterozygosity_pct\trepeat_pct',
                     newLine: true)
        .ifEmpty(file('NO_GENOMESCOPE_METRICS'))


    // ---- Per-sample run info (evidence + strategy) for the report header ----
    ch_run_info_tsv = ch_input
        .map { meta, reads -> tuple(meta.sample, meta) }
        .join(ch_ploidy_by_sample)
        .map { sample, meta, ploidy ->
            [ meta.sample, meta.hifi, meta.hic, meta.tellseq, meta.shortread,
              meta.assembler, meta.n_hap, ploidy, meta.dedup, meta.mito_tool ]
                .collect { it == null ? '' : it }.join('\t')
        }
        .collectFile(name: 'run_info.tsv',
                     seed: ['sample','hifi','hic','tellseq','shortread',
                            'assembler','n_hap','ploidy','dedup','mito_tool'].join('\t'),
                     newLine: true)
        .ifEmpty(file('NO_RUN_INFO'))

    // ---- Workflow provenance for the report header ----
    // Render to plain Strings HERE. workflow.start is a java.time.OffsetDateTime and
    // nextflow.version is a VersionNumber; if either reaches a channel still inside a
    // GString, Kryo can't serialize it ("Unable to create serializer ... OffsetDateTime").
    // .toString() on each line collapses the embedded objects to their string form first.
    def wf_lines = [
            "key\tvalue",
            "pipeline\t${workflow.manifest.name   ?: 'gcl_genome_assembly'}",
            "version\t${workflow.manifest.version ?: 'unknown'}",
            "revision\t${workflow.revision        ?: 'unknown'}",
            "commit\t${workflow.commitId          ?: 'unknown'}",
            "run_name\t${workflow.runName         ?: 'unknown'}",
            "profile\t${workflow.profile          ?: 'unknown'}",
            "nextflow\t${nextflow.version}",
            "start\t${workflow.start}"
        ].collect { it.toString() }.join('\n')

    ch_workflow_info = Channel
        .of(wf_lines)
        .collectFile(name: 'workflow_info.tsv', newLine: true)

    COLLECT_SOFTWARE_VERSIONS(
        ch_versions.collectFile(name: 'software_versions_raw.tsv')
    )

    // ---- Call SUMMARY_REPORT ----
    SUMMARY_REPORT(
        ch_report_manifest,
        ch_compiled_metrics,
        ch_telomere_for_report,
        ch_pairwise_summary,
        ch_mito_stats_for_report,
        ch_teloclip_stats_for_report,
        ch_sample_taxonomy_tsv,
        ch_genome_size_tsv,
        ch_genomescope_metrics_tsv,
        ch_workflow_info,
        ch_run_info_tsv,
        ch_pangenome_report,
        ch_name_map,
        COLLECT_SOFTWARE_VERSIONS.out.versions,
        ch_summary_report_script
    )
}
