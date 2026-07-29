/*
    ANNOTATE_PLASTID — plant plastome annotation (Chloë). Drops sub-<min_len> records,
    array-guards the GFF glob (never grabs the input), degrades to a note on failure.
*/
process ANNOTATE_PLASTID {
    tag "${meta.sample}:${org_type}"
    label 'chloe'
    publishDir "${params.outdir}/organelle/annotation/${org_type}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(org_type), path(fasta)
    val chloe_dir

    output:
    tuple val(meta), val(org_type), path("${meta.sample}.${org_type}.chloe.gff"),           emit: gff,  optional: true
    tuple val(meta), val(org_type), path("${meta.sample}.${org_type}.annotation_note.txt"), emit: note, optional: true
    tuple val(meta), val(org_type), path("chloe_raw"),                                       emit: raw,  optional: true
    path "versions.tsv", emit: versions

    script:
    def min_len = params.chloe_min_contig ?: 200
    """
    set -eu
    shopt -s nullglob
    export JULIA_DEPOT_PATH="${chloe_dir}/depot"
    SAMPLE="${meta.sample}"
    ORG="${org_type}"

    awk '/^>/{h=\$0; next}{s[h]=s[h]\$0} END{for(k in s) if(length(s[k])>=${min_len}) printf "%s\\n%s\\n", k, s[k]}' "${fasta}" > input.fa
    NREC=\$(grep -c '^>' input.fa || echo 0)
    echo "[ANNOTATE_PLASTID] \${NREC} record(s) >= ${min_len} bp fed to Chloe"

    julia --project="${chloe_dir}/chloe" "${chloe_dir}/chloe/chloe.jl" annotate input.fa \\
        || echo "[ANNOTATE_PLASTID] chloe exited non-zero"

    mkdir -p chloe_raw
    produced=( *.chloe.gff *.chloe.sff *.gff *.sff )
    [ \${#produced[@]} -gt 0 ] && cp "\${produced[@]}" chloe_raw/ 2>/dev/null || true

    gffs=( *.chloe.gff *.gff )
    if [ \${#gffs[@]} -gt 0 ]; then
        cp "\$(ls -S "\${gffs[@]}" | head -n1)" "\${SAMPLE}.\${ORG}.chloe.gff"
        echo "[ANNOTATE_PLASTID] annotation written"
    else
        cat > "\${SAMPLE}.\${ORG}.annotation_note.txt" <<'EOF'
Plastome assembled but Chloe produced no annotation. Chloe expects a complete/circular plastome;
a fragmented assembly can crash its ORF finder. Options: improve plastome contiguity, or annotate
manually with GeSeq (web): https://chlorobox.mpimp-golm.mpg.de/geseq.html
EOF
        echo "[ANNOTATE_PLASTID] no GFF produced; wrote annotation_note.txt"
    fi

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
