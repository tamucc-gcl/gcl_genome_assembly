/*
========================================================================================
    ANNOTATE_MITO — metazoan / fungal mitogenome annotation (MITOS2)  ->  MitoHiFi-style output
========================================================================================
    Repo location: modules/annotate_mito.nf

    Runs MITOS2 (runmitos), then converts its BED annotation to a GenBank file, and emits the
    SAME file set MitoHiFi produces so downstream (report, circular map, filtering) treats a
    short-read mitogenome identically to a HiFi one:

      <sample>_<org>_mitogenome.fasta   the assembly (== MitoHiFi mitogenome)
      <sample>_<org>_mitogenome.gb      GenBank from MITOS2 (== MitoHiFi annotation) -> feeds the circular map
      <sample>_<org>_mito_stats.tsv     same columns + counting as MitoHiFi's stats
      <sample>_<org>_mito_contigs.fasta candidate contigs (== the assembly here)
      <sample>_<org>_gene_map.png       MITOS linear gene diagram (optional)
      <sample>_<org>.mitos.bed          raw MITOS annotation (optional)

    Requirements: conda mitos=2.1.10 (binary is `runmitos`, bundles biopython for the converter);
    RefSeq set pre-fetched by DOWNLOAD_MITOS_DB (refseq89m animal / refseq89f fungus).
    MITOS annotates ONE sequence per run — a closed metazoan mito is a single record.

    Input:
      tuple(meta, org_type, status, gcode, refseq, fasta)   status in {circular,linear}
      val  mitos_db     parent dir for the refseq set(s)
      path gb_script    mitos_to_genbank.py
========================================================================================
*/

process ANNOTATE_MITO {
    tag "${meta.sample}:${org_type}"
    label 'mitos'
    publishDir "${params.outdir}/organelle", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(org_type), val(status), val(gcode), val(refseq), path(fasta)
    val mitos_db
    path gb_script

    output:
    tuple val(meta), val(org_type), path("${meta.sample}_${org_type}_mitogenome.gb"),      emit: annotation
    tuple val(meta), val(org_type), path("${meta.sample}_${org_type}_mito_stats.tsv"),     emit: stats
    tuple val(meta), val(org_type), path("${meta.sample}_${org_type}_gene_map.png"),       emit: gene_map, optional: true
    tuple val(meta), val(org_type), path("${meta.sample}_${org_type}.mitos.bed"),          emit: bed,      optional: true
    path "versions.tsv", emit: versions

    script:
    def linear = (status == 'circular') ? '' : '--linear'
    def circ   = (status == 'circular') ? '--circular' : ''
    """
    set -eu
    SAMPLE="${meta.sample}"
    ORG="${org_type}"

    mkdir -p mitos_out
    runmitos \\
        -i "${fasta}" \\
        -c ${gcode} \\
        -o mitos_out \\
        -R "${mitos_db}" \\
        -r ${refseq} \\
        ${linear} \\
        --noplots || echo "[ANNOTATE_MITO] runmitos exited non-zero; collecting whatever landed"

    BED=\$(find mitos_out -name 'result.bed' | head -n1 || true)
    PNG=\$(find mitos_out -name 'result.png' | head -n1 || true)

    [ -n "\${BED}" ] && cp "\${BED}" "\${SAMPLE}_\${ORG}.mitos.bed" || true
    [ -n "\${PNG}" ] && cp "\${PNG}" "\${SAMPLE}_\${ORG}_gene_map.png" || true

    # MITOS BED -> GenBank (single record, CDS/tRNA/rRNA + /gene + /product) for the circular map + report
    GB="\${SAMPLE}_\${ORG}_mitogenome.gb"
    if [ -n "\${BED}" ]; then
        python "${gb_script}" \\
            --fasta "${fasta}" \\
            --bed "\${BED}" \\
            --sample_id "\${SAMPLE}" \\
            --organelle "\${ORG}" \\
            --genetic_code ${gcode} \\
            ${circ} \\
            --output "\${GB}"
    else
        echo "[ANNOTATE_MITO] no result.bed; emitting empty GenBank"
        : > "\${GB}"
    fi

    # stats — SAME columns and counting rules as MitoHiFi's mito_stats.tsv
    MITO_LEN=\$(grep -v '^>' "${fasta}" | tr -d '\\n' | wc -c)
    CIRCULAR=\$([ "${status}" = "circular" ] && echo yes || echo no)
    GENE_COUNT=\$(grep -c '^ */gene=' "\${GB}" 2>/dev/null || true)
    TRNA_COUNT=\$(grep -c '/product="tRNA' "\${GB}" 2>/dev/null || true)
    RRNA_COUNT=\$(grep -c 'rRNA' "\${GB}" 2>/dev/null || true)
    printf 'sample_id\\tmitogenome_length\\tcircular\\tgene_count\\ttrna_count\\trrna_count\\tgenetic_code\\n' > "\${SAMPLE}_\${ORG}_mito_stats.tsv"
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "\${SAMPLE}" "\${MITO_LEN}" "\${CIRCULAR}" "\${GENE_COUNT:-0}" "\${TRNA_COUNT:-0}" "\${RRNA_COUNT:-0}" "${gcode}" >> "\${SAMPLE}_\${ORG}_mito_stats.tsv"
    echo "[ANNOTATE_MITO] \${SAMPLE} \${ORG}: len=\${MITO_LEN} circular=\${CIRCULAR} genes=\${GENE_COUNT:-0} tRNA=\${TRNA_COUNT:-0} rRNA=\${RRNA_COUNT:-0}"

    printf 'MITOS2\\t%s\\n' "\$(runmitos --version 2>&1 | head -n1 || echo NA)" > versions.tsv
    """

    stub:
    """
    printf 'LOCUS       stub\\n//\\n' > ${meta.sample}_${org_type}_mitogenome.gb
    printf 'sample_id\\tmitogenome_length\\tcircular\\tgene_count\\ttrna_count\\trrna_count\\tgenetic_code\\n%s\\t16000\\tyes\\t13\\t22\\t2\\t2\\n' "${meta.sample}" > ${meta.sample}_${org_type}_mito_stats.tsv
    touch versions.tsv
    """
}
