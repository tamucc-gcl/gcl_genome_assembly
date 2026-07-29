/*
========================================================================================
    ORGANELLE_ANNOTATION — annotate assembled organelles (selector subworkflow)
========================================================================================
    Repo location: workflows/organelle_annotation.nf

    Routes each organelle:
      animal_mt / fungus_mt     -> MITOS2 (ANNOTATE_MITO) -> MitoHiFi-style outputs; circular
                                   ones additionally get a circular gene map via MITO_CIRCULAR_MAP
                                   (the same module/script the HiFi path uses)
      embplant_pt / embplant_mt -> GeSeq note (no trusted CLI annotator)

    Emits the same channel shapes as the MitoHiFi path (mitogenome / annotation / stats /
    contigs / gene_map / circular_map) so a short-read mito is interchangeable downstream,
    plus `notes` for the plant organelles.

    take:
      ch_assembly  tuple(meta, org_type, fasta)     SHORTREAD_ORGANELLE.out.assembly
      ch_stats     tuple(meta, org_type, stats_tsv) SHORTREAD_ORGANELLE.out.stats
      ch_gcode     tuple(taxid, genetic_code)
========================================================================================
*/

include { DOWNLOAD_MITOS_DB } from '../modules/download_mitos_db.nf'
include { ANNOTATE_MITO }     from '../modules/annotate_mito.nf'
include { MITO_CIRCULAR_MAP } from '../modules/mito_circular_map.nf'
include { mitosRefseqFor }    from '../functions/taxonomy.nf'

// Scripts (params let the standalone test point at the repo's py_scripts/; default is production layout).
mito_gb_script       = file(params.mitos_to_genbank_script ?: "${projectDir}/py_scripts/mitos_to_genbank.py", checkIfExists: true)
mito_circular_script = file(params.mito_circular_script    ?: "${projectDir}/py_scripts/plot_mito_circular.py", checkIfExists: true)

// Plant organelles: no CLI annotator we trust -> GeSeq note (handles both plant plastid + mito).
process ANNOTATION_NOTE {
    tag "${meta.sample}:${org_type}"
    label 'tiny'
    publishDir "${params.outdir}/organelle/annotation/${org_type}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(org_type)

    output:
    tuple val(meta), val(org_type), path("${meta.sample}.${org_type}.annotation_note.txt"), emit: note

    script:
    def what = (org_type == 'embplant_pt') ? 'Plastome' : 'Plant mitochondrion'
    """
    cat > "${meta.sample}.${org_type}.annotation_note.txt" <<'EOF'
${what} assembled but not automatically annotated: no maintained command-line annotator is
suitable for it in this pipeline. Annotate manually with GeSeq (web-only), which handles both
plant plastid and plant mitochondrial genomes: https://chlorobox.mpimp-golm.mpg.de/geseq.html
EOF
    """

    stub:
    """
    echo "stub note" > ${meta.sample}.${org_type}.annotation_note.txt
    """
}

workflow ORGANELLE_ANNOTATION {
    take:
    ch_assembly
    ch_stats
    ch_gcode

    main:
    // status (circular/linear) per assembly, keyed by (sample, org_type)
    ch_status = ch_stats.map { meta, org, tsv ->
        def cols = tsv.text.trim().readLines().last().split('\t')
        tuple( [meta.sample, org], (cols.size() > 2 ? cols[2] : 'unknown') )
    }

    ch_asm = ch_assembly
        .map { meta, org, fa -> tuple([meta.sample, org], meta, org, fa) }
        .join( ch_status )
        .map { key, meta, org, fa, status -> tuple(meta, org, status, fa) }

    ch_asm
        .branch { meta, org, status, fa ->
            mito: org.endsWith('_mt') && org != 'embplant_mt'          // animal_mt, fungus_mt
            note: org == 'embplant_mt' || org == 'embplant_pt'         // plant organelles -> GeSeq note
        }
        .set { br }

    // -- metazoan / fungal mito -> MITOS2 (MitoHiFi-style outputs) --
    ch_mitos_db = DOWNLOAD_MITOS_DB(
                      br.mito.map { meta, org, status, fa -> params.mitos_downloads }.unique(),
                      params.mitos_refseq_sets,
                      params.mitos_force_download
                  ).db.first()

    br.mito
        .map { meta, org, status, fa -> tuple(meta.taxid?.toString(), meta, org, status, fa) }
        .combine( ch_gcode, by: 0 )
        .map { taxid, meta, org, status, fa, gcode ->
               tuple(meta, org, status, gcode, (params.mitos_refseq ?: mitosRefseqFor(org)), fa) }
        .set { ch_mito_in }
    ANNOTATE_MITO( ch_mito_in, ch_mitos_db, mito_gb_script )

    // circular mitos -> circular gene map (reuse the HiFi module + script). Join status back on.
    ANNOTATE_MITO.out.annotation
        .map { meta, org, gb -> tuple([meta.sample, org], meta, gb) }
        .join( br.mito.map { meta, org, status, fa -> tuple([meta.sample, org], status) } )
        .filter { key, meta, gb, status -> status == 'circular' }
        .map { key, meta, gb, status -> tuple(meta, gb) }
        .set { ch_circular_gb }
    MITO_CIRCULAR_MAP( ch_circular_gb, mito_circular_script )

    // -- plant plastid + plant mito -> GeSeq note --
    ANNOTATION_NOTE( br.note.map { meta, org, status, fa -> tuple(meta, org) } )

    emit:
    mitogenome   = ANNOTATE_MITO.out.mitogenome
    annotation   = ANNOTATE_MITO.out.annotation
    mito_stats   = ANNOTATE_MITO.out.stats
    mito_contigs = ANNOTATE_MITO.out.contigs
    gene_map     = ANNOTATE_MITO.out.gene_map
    circular_map = MITO_CIRCULAR_MAP.out.circular_map
    notes        = ANNOTATION_NOTE.out.note
    versions     = ANNOTATE_MITO.out.versions
}
