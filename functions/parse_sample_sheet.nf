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
========================================================================================
*/

def parseSampleSheet(sample_sheet_path) {
    
    // Parse the sample sheet into a list first
    // This allows each sample to be cached independently
    def samples = []
    def sample_sheet_file = file(sample_sheet_path)
    
    if (!sample_sheet_file.exists()) {
        error "Sample sheet does not exist: ${sample_sheet_path}"
    }
    
    def lines = sample_sheet_file.readLines()
    
    // Process each line
    lines.eachWithIndex { line, idx ->
        // Skip header (assume first line or lines starting with 'sample_id')
        if (idx == 0 || line.startsWith('sample_id')) {
            log.info "Parsing sample sheet: ${sample_sheet_path}"
            return
        }
        
        // Skip empty lines
        if (!line.trim()) return
        
        // Parse CSV line
        def fields = line.split(',').collect { it.trim() }
        
        // Check for correct number of fields
        if (fields.size() != 4) {
            log.warn "Line ${idx + 1}: Malformed CSV line (expected 4 fields, got ${fields.size()}): ${line}"
            return
        }
        
        def (sample_id, hifi_bam_path, hic_r1_path, hic_r2_path) = fields
        
        // Check for missing required fields
        def missing = []
        if (!sample_id) missing << "sample_id"
        if (!hifi_bam_path) missing << "hifi_bam"
        if (!hic_r1_path) missing << "hic_r1"
        if (!hic_r2_path) missing << "hic_r2"
        
        // Skip row if any required fields are missing
        if (missing.size() > 0) {
            log.warn "Line ${idx + 1}: Skipping sample '${sample_id ?: 'UNKNOWN'}': missing required field(s): ${missing.join(', ')}"
            return
        }
        
        // Convert to file objects and check existence
        def hifi_bam = file(hifi_bam_path)
        def hic_r1 = file(hic_r1_path)
        def hic_r2 = file(hic_r2_path)
        
        if (!hifi_bam.exists()) {
            log.warn "Line ${idx + 1}: Skipping sample '${sample_id}': hifi_bam file does not exist: ${hifi_bam_path}"
            return
        }
        if (!hic_r1.exists()) {
            log.warn "Line ${idx + 1}: Skipping sample '${sample_id}': hic_r1 file does not exist: ${hic_r1_path}"
            return
        }
        if (!hic_r2.exists()) {
            log.warn "Line ${idx + 1}: Skipping sample '${sample_id}': hic_r2 file does not exist: ${hic_r2_path}"
            return
        }
        
        // Add validated sample to list
        samples << tuple(sample_id, hifi_bam, hic_r1, hic_r2)
        log.info "Line ${idx + 1}: Successfully parsed sample '${sample_id}'"
    }
    
    // Report results
    if (samples.size() == 0) {
        error "No valid samples found in sample sheet: ${sample_sheet_path}"
    }
    
    log.info "Successfully parsed ${samples.size()} sample(s) from sample sheet"
    
    // CRITICAL: Use Channel.fromList() to create independent cache entries per sample
    // This ensures that adding new samples won't invalidate cache for existing samples
    return Channel.fromList(samples)
}