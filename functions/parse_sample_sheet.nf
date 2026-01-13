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
========================================================================================
*/

def parseSampleSheet(sample_sheet_path) {
    
    return Channel
        .fromPath(sample_sheet_path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            // Validate required columns exist
            if (!row.sample_id) {
                error "Missing 'sample_id' column in sample sheet"
            }
            if (!row.hifi_bam) {
                error "Missing 'hifi_bam' column in sample sheet for sample: ${row.sample_id}"
            }
            if (!row.hic_r1) {
                error "Missing 'hic_r1' column in sample sheet for sample: ${row.sample_id}"
            }
            if (!row.hic_r2) {
                error "Missing 'hic_r2' column in sample sheet for sample: ${row.sample_id}"
            }
            
            // Validate files exist
            def hifi_bam = file(row.hifi_bam, checkIfExists: true)
            def hic_r1 = file(row.hic_r1, checkIfExists: true)
            def hic_r2 = file(row.hic_r2, checkIfExists: true)
            
            // Return tuple with validated data
            tuple(
                row.sample_id,
                hifi_bam,
                hic_r1,
                hic_r2
            )
        }
}