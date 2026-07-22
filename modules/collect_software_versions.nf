process COLLECT_SOFTWARE_VERSIONS {
    tag "software_versions"
    label 'summarize_assembly'
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

    input:
    path 'versions/*'          // per-process versions.tsv files (lines: tool<TAB>version)

    output:
    path 'software_versions.tsv', emit: versions

    script:
    """
    printf 'tool\\tversion\\n' > software_versions.tsv
    cat versions/* 2>/dev/null | sed '/^[[:space:]]*\$/d' | sort -u >> software_versions.tsv
    """
}