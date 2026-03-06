/*
========================================================================================
    PARSE SAMPLE SHEET FUNCTION (CACHE-OPTIMIZED)
========================================================================================
    Reads and validates the input sample sheet CSV file
    
    Expected CSV columns:
        - sample_id: Unique sample identifier
        - hifi_bam: Path to HiFi BAM file
        - hic_r1: Path to Hi-C R1 FASTQ file
        - hic_r2: Path to Hi-C R2 FASTQ file
    
    Rows with missing information are skipped and logged
    
    CACHING BEHAVIOR:
    This function creates independent channel items per sample using .fromList(),
    which ensures that adding new samples to the samplesheet will NOT invalidate
    the cache for existing samples when using -resume.
    
    PATH RESOLUTION:
    Relative paths in the samplesheet are resolved against the samplesheet's
    parent directory. This means paths like "raw_bam/foo.bam" in a samplesheet
    at "data/assembly_samplesheet.csv" resolve to "data/raw_bam/foo.bam".
    Absolute paths (starting with "/") are used as-is.
    
    This approach allows the pipeline to be launched from any directory
    (e.g., .nf/assembly/) as long as the samplesheet path itself is absolute.
========================================================================================
*/

def parseSampleSheet(sample_sheet_path) {

    def samples = []
    def sample_sheet_file = file(sample_sheet_path)

    if (!sample_sheet_file.exists()) {
        error "Sample sheet does not exist: ${sample_sheet_path}"
    }

    // Resolve relative paths against the samplesheet's parent directory
    def sheet_dir = sample_sheet_file.parent
    log.info "Resolving samplesheet paths relative to: ${sheet_dir}"

    def resolvePath = { String p ->
        if (!p) return null
        return p.startsWith('/') ? file(p) : file("${sheet_dir}/${p}")
    }

    def lines = sample_sheet_file.readLines()

    lines.eachWithIndex { line, idx ->
        if (idx == 0 || line.startsWith('sample_id')) {
            log.info "Parsing sample sheet: ${sample_sheet_path}"
            return
        }
        if (!line.trim()) return

        def fields = line.split(',').collect { it.trim() }
        if (fields.size() != 4) {
            log.warn "Line ${idx + 1}: Malformed CSV line (expected 4 fields, got ${fields.size()}): ${line}"
            return
        }

        def (sample_id, hifi_bam_path, hic_r1_path, hic_r2_path) = fields

        def missing = []
        if (!sample_id) missing << "sample_id"
        if (!hifi_bam_path) missing << "hifi_bam"
        if (!hic_r1_path) missing << "hic_r1"
        if (!hic_r2_path) missing << "hic_r2"
        if (missing) {
            log.warn "Line ${idx + 1}: Skipping sample '${sample_id ?: 'UNKNOWN'}': missing required field(s): ${missing.join(', ')}"
            return
        }

        def hifi_bam = resolvePath(hifi_bam_path)
        def hic_r1   = resolvePath(hic_r1_path)
        def hic_r2   = resolvePath(hic_r2_path)

        if (!hifi_bam.exists()) {
            log.warn "Line ${idx + 1}: Skipping sample '${sample_id}': hifi_bam file does not exist: ${hifi_bam_path} (resolved to: ${hifi_bam})"
            return
        }
        if (!hic_r1.exists()) {
            log.warn "Line ${idx + 1}: Skipping sample '${sample_id}': hic_r1 file does not exist: ${hic_r1_path} (resolved to: ${hic_r1})"
            return
        }
        if (!hic_r2.exists()) {
            log.warn "Line ${idx + 1}: Skipping sample '${sample_id}': hic_r2 file does not exist: ${hic_r2_path} (resolved to: ${hic_r2})"
            return
        }

        samples << tuple(sample_id, hifi_bam, hic_r1, hic_r2)
        log.info "Line ${idx + 1}: Successfully parsed sample '${sample_id}'"
    }

    if (samples.size() == 0) {
        error "No valid samples found in sample sheet: ${sample_sheet_path}"
    }

    log.info "Successfully parsed ${samples.size()} sample(s) from sample sheet"
    return Channel.fromList(samples)
}