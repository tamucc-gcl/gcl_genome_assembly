/*
========================================================================================
    GENERATE DECONTAMINATION EVIDENCE WORKFLOW
========================================================================================
    Purpose:
    - Generate supporting evidence for decontamination decisions
    - Coverage evidence: map HiFi reads to cleaned assemblies
    - Taxonomy evidence: DIAMOND blastx against protein database
    - Visualization: BlobTools2 plots and reports
    - Fully parallelized across haplotypes
    
    Design:
    - Takes cleaned assemblies and databases as input
    - Each haplotype generates evidence independently
    - Combines coverage + taxonomy + FCS reports
    - Optional workflow - can be skipped
========================================================================================
*/

nextflow.enable.dsl=2

include { MAP_READS_MINIMAP2 }       from '../modules/map_reads_minimap2.nf'
include { DIAMOND_BLASTX }           from '../modules/diamond_blastx.nf'
include { BLOBTOOLS_CREATE }         from '../modules/blobtools2_create.nf'
include { BLOBTOOLS_VIEWPLOT }       from '../modules/blobtools2_viewplot.nf'
include { FCS_BLOB_EVIDENCE_REPORT } from '../modules/fcs_blob_evidence_report.nf'

workflow GENERATE_DECONTAM_EVIDENCE {
    take:
    decontaminated   // channel: tuple(haplotype_id, clean_fasta)
    contaminants     // channel: tuple(haplotype_id, contam_fasta)
    action_reports   // channel: tuple(haplotype_id, action_report)
    taxonomy_reports // channel: tuple(haplotype_id, taxonomy_report)
    hifi_reads       // channel: tuple(sample_id, hifi_fastq)
    diamond_db       // channel: DIAMOND database path (already a channel)
    taxdump_dir      // channel: NCBI taxonomy directory (already a channel)
    
    main:
    
    /*
    ========================================================================================
        STEP 1: Map HiFi Reads to Cleaned Assemblies (Coverage Evidence)
    ========================================================================================
    */
    
    // Extract sample_id from haplotype_id and join with HiFi reads
    decontaminated
        .map { haplotype_id, clean_fasta ->
            def sample_id = haplotype_id.replaceAll(/_hap[12]$/, '')
            tuple(sample_id, haplotype_id, clean_fasta)
        }
        .combine(hifi_reads, by: 0)
        .set { ch_combined }
    
    // Split into separate channels for each input
    ch_combined
        .map { sample_id, haplotype_id, clean_fasta, hifi_fastq -> clean_fasta }
        .set { ch_assembly }
    
    ch_combined
        .map { sample_id, haplotype_id, clean_fasta, hifi_fastq -> hifi_fastq }
        .set { ch_reads }
    
    ch_combined
        .map { sample_id, haplotype_id, clean_fasta, hifi_fastq -> params.evidence.map_preset ?: 'map-hifi' }
        .set { ch_preset }
    
    MAP_READS_MINIMAP2(ch_assembly, ch_reads, ch_preset)
    
    // Restore haplotype_id for BAM output
    decontaminated
        .map { haplotype_id, clean_fasta -> haplotype_id }
        .combine(MAP_READS_MINIMAP2.out.bam)
        .set { ch_bam_with_id }
    
    /*
    ========================================================================================
        STEP 2: DIAMOND BLASTX (Taxonomy Evidence)
    ========================================================================================
    */
    decontaminated
        .map { haplotype_id, clean_fasta -> tuple(haplotype_id, clean_fasta) }
        .combine(diamond_db)  // diamond_db is already a channel, don't wrap in Channel.value()
        .set { ch_diamond_combined }
    
    ch_diamond_combined
        .map { haplotype_id, clean_fasta, dmnd -> clean_fasta }
        .set { ch_diamond_assembly }
    
    ch_diamond_combined
        .map { haplotype_id, clean_fasta, dmnd -> dmnd }
        .set { ch_diamond_db }
    
    ch_diamond_combined
        .map { haplotype_id, clean_fasta, dmnd -> params.evidence.diamond_max_target_seqs ?: 1 }
        .set { ch_diamond_max_seqs }
    
    ch_diamond_combined
        .map { haplotype_id, clean_fasta, dmnd -> params.evidence.diamond_evalue ?: 1e-25 }
        .set { ch_diamond_evalue }
    
    DIAMOND_BLASTX(ch_diamond_assembly, ch_diamond_db, ch_diamond_max_seqs, ch_diamond_evalue)
    
    // Restore haplotype_id for DIAMOND output
    decontaminated
        .map { haplotype_id, clean_fasta -> haplotype_id }
        .combine(DIAMOND_BLASTX.out.out_hits)
        .set { ch_diamond_with_id }
    
    /*
    ========================================================================================
        STEP 3: BlobTools2 Create + Visualize
    ========================================================================================
    */
    decontaminated
        .map { haplotype_id, clean_fasta -> tuple(haplotype_id, clean_fasta) }
        .join(ch_diamond_with_id, by: 0)
        .join(ch_bam_with_id, by: 0)
        .combine(taxdump_dir)  // taxdump_dir is already a channel, don't wrap in Channel.value()
        .set { ch_blob_combined }
    
    ch_blob_combined
        .map { haplotype_id, clean_fasta, hits, bam, taxdump -> clean_fasta }
        .set { ch_blob_assembly }
    
    ch_blob_combined
        .map { haplotype_id, clean_fasta, hits, bam, taxdump -> hits }
        .set { ch_blob_hits }
    
    ch_blob_combined
        .map { haplotype_id, clean_fasta, hits, bam, taxdump -> bam }
        .set { ch_blob_bam }
    
    ch_blob_combined
        .map { haplotype_id, clean_fasta, hits, bam, taxdump -> taxdump }
        .set { ch_blob_taxdump }
    
    ch_blob_combined
        .map { haplotype_id, clean_fasta, hits, bam, taxdump -> params.evidence.blob_min_contig_len ?: 1000 }
        .set { ch_blob_min_len }
    
    BLOBTOOLS_CREATE(ch_blob_assembly, ch_blob_hits, ch_blob_bam, ch_blob_taxdump, ch_blob_min_len)
    
    // Restore haplotype_id for blobtools output
    decontaminated
        .map { haplotype_id, clean_fasta -> haplotype_id }
        .combine(BLOBTOOLS_CREATE.out.out_blobdir)
        .set { ch_blobdir_with_id }
    
    BLOBTOOLS_VIEWPLOT(BLOBTOOLS_CREATE.out.out_blobdir)
    
    // Restore haplotype_id for plots output
    decontaminated
        .map { haplotype_id, clean_fasta -> haplotype_id }
        .combine(BLOBTOOLS_VIEWPLOT.out.out_dir)
        .set { ch_plots_with_id }
    
    /*
    ========================================================================================
        STEP 4: Generate Comprehensive Evidence Report
    ========================================================================================
    */
    action_reports
        .join(taxonomy_reports, by: 0)
        .join(ch_plots_with_id, by: 0)
        .join(decontaminated, by: 0)
        .join(contaminants, by: 0)
        .set { ch_report_combined }
    
    ch_report_combined
        .map { haplotype_id, action_report, taxonomy_report, plots_dir, clean_fasta, contam_fasta -> action_report }
        .set { ch_report_action }
    
    ch_report_combined
        .map { haplotype_id, action_report, taxonomy_report, plots_dir, clean_fasta, contam_fasta -> taxonomy_report }
        .set { ch_report_taxonomy }
    
    ch_report_combined
        .map { haplotype_id, action_report, taxonomy_report, plots_dir, clean_fasta, contam_fasta -> plots_dir }
        .set { ch_report_plots }
    
    ch_report_combined
        .map { haplotype_id, action_report, taxonomy_report, plots_dir, clean_fasta, contam_fasta -> clean_fasta }
        .set { ch_report_clean }
    
    ch_report_combined
        .map { haplotype_id, action_report, taxonomy_report, plots_dir, clean_fasta, contam_fasta -> contam_fasta }
        .set { ch_report_contam }
    
    FCS_BLOB_EVIDENCE_REPORT(ch_report_action, ch_report_taxonomy, ch_report_plots, ch_report_clean, ch_report_contam)
    
    // Restore haplotype_id for report outputs
    decontaminated
        .map { haplotype_id, clean_fasta -> haplotype_id }
        .combine(FCS_BLOB_EVIDENCE_REPORT.out.actions_tsv)
        .set { ch_actions_with_id }
    
    decontaminated
        .map { haplotype_id, clean_fasta -> haplotype_id }
        .combine(FCS_BLOB_EVIDENCE_REPORT.out.annotated_tsv)
        .set { ch_annotated_with_id }
    
    decontaminated
        .map { haplotype_id, clean_fasta -> haplotype_id }
        .combine(FCS_BLOB_EVIDENCE_REPORT.out.report_md)
        .set { ch_report_with_id }
    
    emit:
    // Coverage evidence
    coverage_bam = ch_bam_with_id             // tuple(haplotype_id, bam)
    
    // Taxonomy evidence
    diamond_hits = ch_diamond_with_id         // tuple(haplotype_id, hits)
    
    // BlobTools outputs
    blobtools_dir = ch_blobdir_with_id        // tuple(haplotype_id, blobdir)
    blobtools_plots = ch_plots_with_id        // tuple(haplotype_id, plots_dir)
    
    // Evidence reports
    actions_table = ch_actions_with_id        // tuple(haplotype_id, actions_tsv)
    annotated_table = ch_annotated_with_id    // tuple(haplotype_id, annotated_tsv)
    evidence_report = ch_report_with_id       // tuple(haplotype_id, report_md)
}