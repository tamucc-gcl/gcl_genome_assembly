/*
========================================================================================
    ANNOTATE_MITO — metazoan / fungal mitogenome annotation (MITOS2)
========================================================================================
    Repo location: modules/annotate_mito.nf

    MITOS2 (runmitos.py) is the standard de-novo annotator for METAZOAN mitogenomes (and
    fungal, via refseq89f). It is NOT appropriate for plant mitochondria — those route
    elsewhere (currently deferred).

    Requirements:
      - conda: mitos>=2  (bioconda)
      - a RefSeq reference set pre-installed from Zenodo (record 2672835 -> newer version),
        e.g. refseq89m (Metazoa) / refseq89f (Fungi), living under params.mitos_downloads.
        Set up once, like the GetOrganelle DB.

    MITOS2 annotates ONE sequence per run. A closed metazoan mitogenome is a single record,
    so we pass the assembly as-is. (A fragmented/multi-record mito would need the longest
    record selected first — a follow-up if one ever shows up here.)

    Invocation is the minimal, cross-version-verified form:
        runmitos.py -i <fasta> -c <code> -o <dir> -R <db_parent> -r <refseqNN> [--linear]
    (No --noplots: the plural/singular of that flag varies by version, so we leave it off
     rather than guess; plots are harmless extra output in the raw dir.)

    Input:
      tuple(meta, org_type, status, gcode, refseq, fasta)
        status ∈ {circular, linear}   -> --linear when not circular
        gcode                          NCBI mito genetic code (geneticCodeFor)
        refseq                         reference set name (refseq89m animal / refseq89f fungus)
      val mitos_db                     parent dir holding the refseq set(s)
    Output:
      bed  tuple(meta, org_type, <sample>.<org>.mitos.bed)   MITOS annotation (optional: absent on failure)
      gff  optional GFF
      faa  optional protein predictions
      raw  full MITOS output dir
========================================================================================
*/

process ANNOTATE_MITO {
    tag "${meta.sample}:${org_type}"
    label 'mitos'
    publishDir "${params.outdir}/organelle/annotation/${org_type}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(org_type), val(status), val(gcode), val(refseq), path(fasta)
    val mitos_db

    output:
    tuple val(meta), val(org_type), path("${meta.sample}.${org_type}.mitos.bed"), emit: bed, optional: true
    tuple val(meta), val(org_type), path("${meta.sample}.${org_type}.mitos.gff"), emit: gff, optional: true
    tuple val(meta), val(org_type), path("${meta.sample}.${org_type}.mitos.faa"), emit: faa, optional: true
    tuple val(meta), val(org_type), path("mitos_out"),                            emit: raw
    path "versions.tsv", emit: versions

    script:
    def linear = (status == 'circular') ? '' : '--linear'
    """
    set -eu
    SAMPLE="${meta.sample}"
    ORG="${org_type}"

    mkdir -p mitos_out
    runmitos.py \\
        -i "${fasta}" \\
        -c ${gcode} \\
        -o mitos_out \\
        -R "${mitos_db}" \\
        -r ${refseq} \\
        ${linear} || echo "[ANNOTATE_MITO] runmitos.py exited non-zero; collecting whatever landed"

    # MITOS writes result.* into -o (sometimes a subdir); collect wherever they land
    bed=\$(find mitos_out -name 'result.bed' | head -n1 || true)
    gff=\$(find mitos_out -name 'result.gff' | head -n1 || true)
    faa=\$(find mitos_out -name 'result.faa' | head -n1 || true)
    [ -n "\${bed}" ] && cp "\${bed}" "\${SAMPLE}.\${ORG}.mitos.bed" || true
    [ -n "\${gff}" ] && cp "\${gff}" "\${SAMPLE}.\${ORG}.mitos.gff" || true
    [ -n "\${faa}" ] && cp "\${faa}" "\${SAMPLE}.\${ORG}.mitos.faa" || true

    printf 'MITOS2\\t%s\\n' "\$(python -c 'import mitos; print(mitos.__version__)' 2>/dev/null || echo NA)" > versions.tsv
    """

    stub:
    """
    mkdir -p mitos_out
    printf '# stub\\n' > ${meta.sample}.${org_type}.mitos.bed
    touch versions.tsv
    """
}
