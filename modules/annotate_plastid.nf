/*
========================================================================================
    ANNOTATE_PLASTID — plant plastome annotation (Chloë)
========================================================================================
    Repo location: modules/annotate_plastid.nf

    Chloë (ian-small/Chloe.jl) is optimised for angiosperm chloroplast genomes. Verified CLI:
        julia --project=<chloe> <chloe>/chloe.jl annotate <fasta>
    Default output is GFF, written next to the input as <input>.chloe.gff. References are
    auto-found as a sibling of the chloe folder (SETUP_CHLOE clones them there).

    We copy the input to a real local file first (cp -L) so Chloë writes its output into the
    task workdir even if Nextflow staged the input as a symlink.

    NOTE: annotation quality tracks assembly contiguity — a fragmented (multi-record) plastome
    will annotate but genes spanning fragment boundaries may come back partial.

    Input:
      tuple(meta, org_type, fasta)
      val chloe_dir   the SETUP_CHLOE install dir
    Output:
      gff  tuple(meta, org_type, <sample>.<org>.chloe.gff)   (optional: absent on failure)
      raw  all Chloë outputs
========================================================================================
*/

process ANNOTATE_PLASTID {
    tag "${meta.sample}:${org_type}"
    label 'chloe'
    publishDir "${params.outdir}/organelle/annotation/${org_type}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(org_type), path(fasta)
    val chloe_dir

    output:
    tuple val(meta), val(org_type), path("${meta.sample}.${org_type}.chloe.gff"), emit: gff, optional: true
    tuple val(meta), val(org_type), path("chloe_raw"),                            emit: raw, optional: true
    path "versions.tsv", emit: versions

    script:
    """
    set -eu
    shopt -s nullglob
    export JULIA_DEPOT_PATH="${chloe_dir}/depot"
    SAMPLE="${meta.sample}"
    ORG="${org_type}"

    # real local copy so Chloë's "write next to input" lands in the workdir
    cp -L "${fasta}" input.fa

    julia --project="${chloe_dir}/chloe" "${chloe_dir}/chloe/chloe.jl" annotate input.fa \\
        || echo "[ANNOTATE_PLASTID] chloe exited non-zero; collecting whatever landed"

    mkdir -p chloe_raw
    produced=( *.chloe.gff *.chloe.sff *.gff *.sff )
    [ \${#produced[@]} -gt 0 ] && cp "\${produced[@]}" chloe_raw/ 2>/dev/null || true

    # canonical annotation = the GFF (prefer *.chloe.gff)
    gff=\$(ls -S *.chloe.gff *.gff 2>/dev/null | head -n1 || true)
    [ -n "\${gff}" ] && cp "\${gff}" "\${SAMPLE}.\${ORG}.chloe.gff" || true

    V=\$(git -C "${chloe_dir}/chloe" rev-parse --short HEAD 2>/dev/null || echo NA)
    printf 'Chloe\\t%s\\n' "\$V" > versions.tsv
    """

    stub:
    """
    printf '##gff-version 3\\n' > ${meta.sample}.${org_type}.chloe.gff
    mkdir -p chloe_raw
    touch versions.tsv
    """
}
