include { MITOHIFI }             from '../modules/mitohifi.nf'
include { MITO_CIRCULAR_MAP }    from '../modules/mito_circular_map.nf'
include { SHORTREAD_ORGANELLE }  from './shortread_organelle.nf'
include { ORGANELLE_ANNOTATION } from './organelle_annotation.nf'

mito_circular_script = file("${projectDir}/py_scripts/plot_mito_circular.py", checkIfExists: true)

workflow ORGANELLE {
    take:
    ch_reads
    ch_mito_ref
    ch_gcode
    ch_organelle
    capabilities

    main:
    ch_assemblies = Channel.empty()
    ch_annotation = Channel.empty()
    ch_stats = Channel.empty()
    ch_circular = Channel.empty()
    ch_gene_map = Channel.empty()
    ch_baits = Channel.empty()
    ch_notes = Channel.empty()
    ch_versions = Channel.empty()

    ch_reads.branch { meta, hifi_fastq, sr1, sr2 ->
        hifi: meta.hifi
        other: true
    }.set { ch_org }

    if (capabilities.hifi) {
        ch_mitohifi_input = ch_org.hifi
            .map { meta, hifi_fastq, sr1, sr2 -> tuple(meta.taxid?.toString(), meta, hifi_fastq) }
            .combine(ch_mito_ref, by: 0)
            .combine(ch_gcode, by: 0)
            .map { taxid, meta, hifi_fastq, ref_fa, ref_gb, gcode -> tuple(meta, hifi_fastq, ref_fa, ref_gb, gcode) }
        MITOHIFI(ch_mitohifi_input)
        MITO_CIRCULAR_MAP(MITOHIFI.out.annotation, mito_circular_script)
        ch_assemblies = MITOHIFI.out.mitogenome
        ch_annotation = MITOHIFI.out.annotation
        ch_stats = MITOHIFI.out.stats
        ch_circular = MITO_CIRCULAR_MAP.out.circular_map
        ch_gene_map = MITOHIFI.out.gene_map
        ch_versions = ch_versions.mix(MITOHIFI.out.versions)
        // Explicit successful completion even when no complete mitogenome was recovered.
        ch_baits = MITOHIFI.out.status
            .map { meta, status -> tuple(meta.sample, status) }
            .join(MITOHIFI.out.mitogenome.map { meta, fa -> tuple(meta.sample, fa) }, remainder: true)
            .map { sample, status, fa -> tuple(sample, fa ? [fa] : []) }
    }

    if (capabilities.shortread_organelle) {
        ch_sr_reads = ch_org.other.map { meta, hifi_fastq, sr1, sr2 -> tuple(meta, sr1, sr2) }
        SHORTREAD_ORGANELLE(ch_sr_reads, ch_organelle)
        ORGANELLE_ANNOTATION(SHORTREAD_ORGANELLE.out.assembly, SHORTREAD_ORGANELLE.out.stats, ch_gcode)
        ch_assemblies = ch_assemblies.mix(SHORTREAD_ORGANELLE.out.assembly.map { meta, org, fa -> tuple(meta, fa) })
        ch_annotation = ch_annotation.mix(ORGANELLE_ANNOTATION.out.annotation.map { meta, org, gb -> tuple(meta, gb) })
        ch_stats = ch_stats.mix(ORGANELLE_ANNOTATION.out.mito_stats.map { meta, org, tsv -> tuple(meta, tsv) })
        ch_circular = ch_circular.mix(ORGANELLE_ANNOTATION.out.circular_map)
        ch_gene_map = ch_gene_map.mix(ORGANELLE_ANNOTATION.out.gene_map.map { meta, org, png -> tuple(meta, png) })
        ch_notes = ORGANELLE_ANNOTATION.out.notes
        ch_versions = ch_versions.mix(SHORTREAD_ORGANELLE.out.versions, ORGANELLE_ANNOTATION.out.versions)
        ch_sr_baits = ch_org.other
            .map { meta, hifi_fastq, sr1, sr2 -> tuple(meta.sample, meta.sample) }
            .join(SHORTREAD_ORGANELLE.out.assembly
                .map { meta, org, fa -> tuple(meta.sample, fa) }.groupTuple(), remainder: true)
            .map { sample, s, fas ->
                tuple(sample, fas == null ? [] : (fas.size() > 1 ? fas.sort { it.name } : fas))
            }
        ch_baits = ch_baits.mix(ch_sr_baits)
    }

    emit:
    assemblies = ch_assemblies
    baits = ch_baits
    annotation = ch_annotation
    stats = ch_stats
    circular_map = ch_circular
    gene_map = ch_gene_map
    notes = ch_notes
    versions = ch_versions
}
