include { MAP_HIFI_FINAL } from '../modules/map_hifi_final.nf'
include { COVERAGE_BOOK } from '../modules/coverage_book.nf'

workflow FINAL_READ_EVIDENCE {
    take:
    assemblies
    hifi
    books_on
    coverage_script

    main:
    assemblies.map { meta, fa -> tuple(meta.sample, meta, fa) }
        .combine(hifi.map { meta, reads -> tuple(meta.sample, reads) }, by: 0)
        .map { sample, meta, fa, reads -> tuple(meta, fa, reads) }
        .set { mapping_input }
    MAP_HIFI_FINAL(mapping_input)
    books = Channel.empty()
    if (books_on) {
        COVERAGE_BOOK(MAP_HIFI_FINAL.out.bam.join(MAP_HIFI_FINAL.out.fai), coverage_script)
        books = COVERAGE_BOOK.out.pdf
    }

    emit:
    bam = MAP_HIFI_FINAL.out.bam
    reference_digest = MAP_HIFI_FINAL.out.reference_digest
    pdf = books
    versions = MAP_HIFI_FINAL.out.versions
}
