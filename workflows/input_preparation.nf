include { INPUT_SUMMARY } from '../modules/input_summary.nf'

workflow INPUT_PREPARATION {
    take:
    parsed

    main:
    ch_samples = Channel.fromList(parsed.samples)
    ch_status = Channel.value(parsed.status)
    INPUT_SUMMARY(ch_status)

    emit:
    samples = ch_samples
    status = ch_status
    report = INPUT_SUMMARY.out.report
}
