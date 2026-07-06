/*
========================================================================================
    PARSE SAMPLE SHEET — wide, header-driven; emits tuple(meta, reads)
========================================================================================
    Repo location: functions/parse_sample_sheet.nf

    Header-driven CSV. A read type is PRESENT for a sample iff its column exists AND the
    cell is non-empty. A missing column, or an empty cell, both mean "no data of that type".

    Recognized columns:
      sample_id                                   (required)
      hifi_bam                                    HiFi reads (BAM)
      hic_r1, hic_r2                              Hi-C paired FASTQ
      tellseq_r1, tellseq_r2                      TellSeq linked-read paired FASTQ
      sr_r1, sr_r2                                short-read shotgun paired FASTQ
      ploidy        (optional)                    'haploid'/'diploid' or a positive integer (organism ploidy)
      n_hap         (optional)                    1 | 2 — OUTPUT haplotype count override
                                                  (blank = derived: spades->1, hifiasm-> ploidy==1?1:2)
      assembler     (optional)                    'hifiasm' | 'spades'
      dedup         (optional)                    'purge_dups' | 'redundans' | 'none'
      mito_tool     (optional)                    'mitohifi' | 'mitofinder' | 'none'
      species       (optional)                    organism scientific name (per-sample; falls back
                                                  to params.mitohifi_species) — organelle/kingdom (4b)
      taxid         (optional)                    NCBI taxonomy id — alternative to species (4b)

    Unrecognized columns are ignored. Optional strategy columns override the matching
    global params (precedence is defined in buildMeta). The legacy 4-column sheet
    (sample_id,hifi_bam,hic_r1,hic_r2) remains valid and resolves to current behavior.

    EMITS: Channel.fromList of tuple(meta, reads)
      meta  : map from buildMeta (functions/meta.nf)
      reads : [ hifi_bam, hic_r1, hic_r2, tellseq_r1, tellseq_r2, sr_r1, sr_r2 ] -> Path | null

    CACHING: uses Channel.fromList — independent per-sample cache entries, so adding a
    sample does not invalidate existing samples under -resume. Paths inside `reads` are
    carried as data; extract the needed path(s) into a top-level `path` input right before
    each process (Nextflow only stages top-level path inputs, not paths nested in a map).

    PATHS: relative paths resolve against the sample sheet's parent directory; absolute as-is.

    ERROR POLICY: per-row problems (bad enum, incomplete pair, missing file, no contig reads)
    log a warning and skip that row (matches prior behavior). Hard error only if no valid
    samples remain, or the header lacks 'sample_id'.
========================================================================================
*/

include { buildMeta } from './meta.nf'

def parseSampleSheet(sample_sheet_path) {

    def sheet = file(sample_sheet_path)
    if (!sheet.exists()) error "Sample sheet does not exist: ${sample_sheet_path}"

    def sheet_dir = sheet.parent
    log.info "Parsing sample sheet: ${sample_sheet_path}"
    log.info "Resolving relative paths against: ${sheet_dir}"

    def resolve = { String p -> (!p) ? null : (p.startsWith('/') ? file(p) : file("${sheet_dir}/${p}")) }

    def lines = sheet.readLines()

    // first non-blank line is the header
    def headerIdx = lines.findIndexOf { it?.trim() }
    if (headerIdx < 0) error "Sample sheet is empty: ${sample_sheet_path}"

    def cols = lines[headerIdx].split(',', -1).collect { it.trim() }
    if (!cols.contains('sample_id'))
        error "Sample sheet header must contain a 'sample_id' column. Found: ${cols.join(', ')}"

    def READ_COLS = ['hifi_bam','hic_r1','hic_r2','tellseq_r1','tellseq_r2','sr_r1','sr_r2']

    def samples = []

    lines.eachWithIndex { line, i ->
        if (i <= headerIdx) return
        if (!line?.trim()) return

        def vals = line.split(',', -1).collect { it.trim() }   // -1 keeps trailing empty cells
        def row = [:]
        cols.eachWithIndex { c, ci -> row[c] = (ci < vals.size()) ? vals[ci] : '' }

        def present = { String c -> (cols.contains(c) && row[c]?.trim()) ? true : false }
        def sid = row['sample_id']

        try {
            def hasHic  = present('hic_r1')     || present('hic_r2')
            def hasTell = present('tellseq_r1') || present('tellseq_r2')
            def hasSr   = present('sr_r1')      || present('sr_r2')

            // paired-end completeness
            if (present('hic_r1')     != present('hic_r2'))     throw new IllegalArgumentException("sample '${sid}': Hi-C needs both hic_r1 and hic_r2")
            if (present('tellseq_r1') != present('tellseq_r2')) throw new IllegalArgumentException("sample '${sid}': TellSeq needs both tellseq_r1 and tellseq_r2")
            if (present('sr_r1')      != present('sr_r2'))      throw new IllegalArgumentException("sample '${sid}': short-read needs both sr_r1 and sr_r2")

            def meta = buildMeta([
                sample_id : sid,
                hifi      : present('hifi_bam'),
                hic       : hasHic,
                tellseq   : hasTell,
                shortread : hasSr,
                assembler : row['assembler'],
                ploidy    : row['ploidy'],
                n_hap     : row['n_hap'],
                dedup     : row['dedup'],
                mito_tool : row['mito_tool'],
                species   : row['species'],
                taxid     : row['taxid']
            ])

            // resolve + existence-check present paths only
            def reads = [:]
            READ_COLS.each { c -> reads[c] = present(c) ? resolve(row[c]) : null }
            def missing = reads.findAll { k, v -> v != null && !v.exists() }
            if (missing)
                throw new IllegalArgumentException("sample '${sid}': file(s) not found: " +
                    missing.collect { k, v -> "${k}=${v}" }.join(', '))

            def readTypes = []
            if (meta.hifi)      readTypes << 'hifi'
            if (meta.hic)       readTypes << 'hic'
            if (meta.tellseq)   readTypes << 'tellseq'
            if (meta.shortread) readTypes << 'shortread'

            samples << tuple(meta, reads)
            log.info "  parsed '${sid}': assembler=${meta.assembler}, " +
                     "ploidy=${meta.ploidy}n -> n_hap=${meta.n_hap}, reads=[${readTypes.join(',')}], " +
                     "dedup=${meta.dedup}, mito=${meta.mito_tool}, species=${meta.species}"
        }
        catch (Exception e) {
            log.warn("Skipping row ${i + 1} (sample '${sid ?: 'UNKNOWN'}'): ${e.class.simpleName}: ${e.message}", e)
        }
    }

    if (samples.size() == 0) error "No valid samples found in sample sheet: ${sample_sheet_path}"
    log.info "Successfully parsed ${samples.size()} sample(s) from ${sample_sheet_path}"
    return Channel.fromList(samples)
}
