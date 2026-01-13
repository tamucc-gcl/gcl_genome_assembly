/*
========================================================================================
    PARSE SAMPLE SHEET FUNCTION
========================================================================================
    Reads and validates the input sample sheet CSV file
    
    Expected CSV columns:
        - sample_id: Unique sample identifier
        - hifi_bam: Path to HiFi BAM file
        - hic_r1: Path to Hi-C R1 FASTQ file
        - hic_r2: Path to Hi-C R2 FASTQ file
    
    Rows with missing information are skipped and logged
========================================================================================
*/

def parseSampleSheet(sample_sheet_path) {
    
    return Channel
        .fromPath(sample_sheet_path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            // Check for missing required columns
            def missing = []
            if (!row.sample_id) missing << "sample_id"
            if (!row.hifi_bam) missing << "hifi_bam"
            if (!row.hic_r1) missing << "hic_r1"
            if (!row.hic_r2) missing << "hic_r2"
            
            // Skip row if any required fields are missing
            if (missing.size() > 0) {
                log.warn "Skipping sample '${row.sample_id ?: 'UNKNOWN'}': missing required field(s): ${missing.join(', ')}"
                return null
            }
            
            // Check if files exist
            def hifi_bam = file(row.hifi_bam)
            def hic_r1 = file(row.hic_r1)
            def hic_r2 = file(row.hic_r2)
            
            if (!hifi_bam.exists()) {
                log.warn "Skipping sample '${row.sample_id}': hifi_bam file does not exist: ${row.hifi_bam}"
                return null
            }
            if (!hic_r1.exists()) {
                log.warn "Skipping sample '${row.sample_id}': hic_r1 file does not exist: ${row.hic_r1}"
                return null
            }
            if (!hic_r2.exists()) {
                log.warn "Skipping sample '${row.sample_id}': hic_r2 file does not exist: ${row.hic_r2}"
                return null
            }
            
            // Return tuple with validated data
            tuple(
                row.sample_id,
                hifi_bam,
                hic_r1,
                hic_r2
            )
        }
        .filter { it != null }  // Remove null entries (skipped samples)
}