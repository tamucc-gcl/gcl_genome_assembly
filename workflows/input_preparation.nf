include { parseSampleSheet } from '../functions/input_validation.nf'
include { INPUT_SUMMARY } from '../modules/input_summary.nf'

workflow INPUT_PREPARATION {
    take:
    sample_sheet
    hic_readsets

    main:
    parsed = parseSampleSheet(sample_sheet, hic_readsets)
    INPUT_SUMMARY(parsed.status.first())

    emit:
    samples = parsed.samples
    status = parsed.status.first()
    report = INPUT_SUMMARY.out.report
}
