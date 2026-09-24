/*
 * READ_PREPARATION owns read conversion, trimming, sample read bundles and
 * genome-size estimation, retaining separate Hi-C library/read-set identities.
 * Input: (meta, reads) per sample; (sample, ploidy) side channel.
 * Moving these processes changes qualified task names and may invalidate cache.
 */
include { BAM_TO_FASTQ } from '../modules/bam_to_fastq.nf'
include { TRIM_HIC } from '../modules/trim_hic.nf'
include { TRIM_SHORTREAD } from '../modules/trim_shortread.nf'
include { ESTIMATE_GENOME_SIZE } from '../modules/estimate_genome_size.nf'

workflow READ_PREPARATION {
    take:
    ch_input
    ch_ploidy_by_sample

    capabilities

    main:
    ch_hifi_fastq = Channel.empty()
    ch_hic_trimmed = Channel.empty()
    ch_versions = Channel.empty()

    if (capabilities.hifi) {
    BAM_TO_FASTQ(
        ch_input.filter { meta, reads -> meta.hifi }
                .map { meta, reads -> tuple(meta, reads.hifi_bam) }
    )
    ch_hifi_fastq = BAM_TO_FASTQ.out.fastq
    }

    
    /*
    ========================================================================================
        STEP 2: Trim Hi-C Reads (extract the Hi-C pair from the reads map)
    ========================================================================================
    */
    ch_hic_raw = ch_input.flatMap { meta, reads ->
        (reads.hic_sets ?: []).collect { rs ->
            tuple([id: "${meta.sample}__${rs.readset_id}", sample: meta.sample,
                   library_id: rs.library_id, readset_id: rs.readset_id,
                   readset_count: reads.hic_sets.size()], rs.r1, rs.r2)
        }
    }
    if (capabilities.hic) {
    TRIM_HIC(ch_hic_raw)
    ch_hic_trimmed = TRIM_HIC.out.trimmed_reads
    }
    // Release a sample as soon as all of ITS read sets finish, in stable ID order.
    ch_hic_prepared = ch_hic_trimmed
        .map { rs, r1, r2 -> tuple(groupKey(rs.sample, rs.readset_count), rs.readset_id, r1, r2) }
        .groupTuple()
        .map { key, ids, r1s, r2s ->
            def order = (0..<ids.size()).toList().sort { ids[it] }
            tuple(key.getGroupTarget(), order.collect { r1s[it] }, order.collect { r2s[it] })
        }
        .join(ch_input.filter { meta, reads -> meta.hic }.map { meta, reads -> tuple(meta.sample, meta) })
        .map { sample, r1s, r2s, meta -> tuple(meta, r1s, r2s) }

    // Optional short-read trimming (fastp): raw shotgun -> adapter/quality-trimmed, or
    // pass-through when off. Feeds the assembly + assembly-QC path; SHORTREAD_QC stays on raw.
    ch_shortread_raw = ch_input
        .filter { meta, reads -> meta.shortread }
        .map    { meta, reads -> tuple(meta, reads.sr_r1, reads.sr_r2) }

    if (capabilities.shortread && params.run_shortread_trim) {
        TRIM_SHORTREAD(ch_shortread_raw)
        ch_shortread_reads = TRIM_SHORTREAD.out.trimmed_reads
        ch_versions = ch_versions.mix(TRIM_SHORTREAD.out.versions)
    } else {
        ch_shortread_reads = ch_shortread_raw
    }

    /*
    ========================================================================================
        STEP 3: Combine HiFi FASTQ with trimmed Hi-C reads
    ========================================================================================
    */
    // Full per-sample read bundle for the selector + organelle + genome-size.
    // remainder:true left-joins keep samples lacking HiFi or Hi-C (null slots), so
    // short-read-only rows flow instead of being dropped by an inner join.
    // Per-modality slots: every sample gets exactly one entry per slot — its processed
    // reads if it has that modality, else a null placeholder from ch_input (immediate).
    // Plain 1:1 joins then emit each sample as soon as ITS OWN reads are ready — no
    // waiting on other samples' BAM_TO_FASTQ / TRIM_HIC / TRIM_SHORTREAD to finish.
    ch_hifi_slot = ch_hifi_fastq
        .map { meta, fq -> [ meta.sample, fq ] }
        .mix( ch_input.filter { meta, reads -> !meta.hifi }.map { meta, reads -> [ meta.sample, null ] } )

    ch_hic_slot = ch_hic_prepared
        .map { meta, r1, r2 -> [ meta.sample, [r1, r2] ] }
        .mix( ch_input.filter { meta, reads -> !meta.hic }.map { meta, reads -> [ meta.sample, null ] } )

    ch_sr_slot = ch_shortread_reads
        .map { meta, r1, r2 -> [ meta.sample, [r1, r2] ] }
        .mix( ch_input.filter { meta, reads -> !meta.shortread }.map { meta, reads -> [ meta.sample, null ] } )

    ch_input
        .map { meta, reads -> [ meta.sample, meta ] }
        .join( ch_hifi_slot )
        .join( ch_hic_slot )
        .join( ch_sr_slot )
        .map { sample, meta, hifi_fastq, hic_pair, sr_pair ->
            def hic_r1 = hic_pair ? hic_pair[0] : null
            def hic_r2 = hic_pair ? hic_pair[1] : null
            def sr_r1  = sr_pair  ? sr_pair[0]  : null
            def sr_r2  = sr_pair  ? sr_pair[1]  : null
            tuple(meta, hifi_fastq, hic_r1, hic_r2, sr_r1, sr_r2)
        }
        .set { ch_reads_all }
        
    // Per-sample reads for assembly QC (meryl DB + mapping): HiFi FASTQ for HiFi samples,
    // the Illumina R1+R2 pair for short-read samples. Read-source-aware QC.
    ch_reads_all
        .map { meta, hifi_fastq, hic_r1, hic_r2, sr_r1, sr_r2 ->
            meta.assembler == 'hifiasm' ? tuple(meta, hifi_fastq) : tuple(meta, [sr_r1, sr_r2])
        }
        .set { ch_qc_reads }

    // Genome-size estimation (jellyfish -> GenomeScope2), concurrent with assembly.
    // Reads by assembler: HiFi for the long-read path, PE for short-read.
    ch_reads_all
        .map { meta, hifi_fastq, r1, r2, sr1, sr2 ->
            def gs_reads = (meta.assembler == 'spades') ? [ sr1, sr2 ] : [ hifi_fastq ]
            tuple(meta.sample, meta, gs_reads)
        }
        .join(ch_ploidy_by_sample)
        .map { sample, meta, gs_reads, ploidy -> tuple(meta, gs_reads, ploidy) }
        .set { ch_gsize_input }

    ESTIMATE_GENOME_SIZE(ch_gsize_input)
    ch_versions = ch_versions.mix(ESTIMATE_GENOME_SIZE.out.versions)

    emit:
    reads = ch_reads_all                       // (meta, hifi, hic_r1, hic_r2, sr_r1, sr_r2)
    qc_reads = ch_qc_reads                     // (meta, read file or paired read list)
    hifi = ch_hifi_fastq               // (meta, fastq)
    hic = ch_hic_prepared                      // (sample meta, ordered R1 list, ordered R2 list)
    hic_raw = ch_hic_raw                       // per read set, for raw QC
    hic_trimmed = ch_hic_trimmed    // per read set, for trimmed QC
    shortread = ch_shortread_reads             // (meta, r1, r2), trimmed or raw
    genome_size = ESTIMATE_GENOME_SIZE.out.size // (meta, size file)
    genome_results = ESTIMATE_GENOME_SIZE.out.results // (meta, summary file)
    versions = ch_versions
}
