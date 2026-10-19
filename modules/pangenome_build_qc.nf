/* Post-build accounting/provenance is separate from the cached Cactus task. */
process PANGENOME_BUILD_QC {
    tag "${taxid}"
    label 'cactus_tools'
    publishDir "${params.outdir}/pangenome/${taxid}/representation", mode: params.publish_dir_mode,
        pattern: 'representation/*', saveAs: { it.replaceFirst('^representation/', '') }
    input:
    tuple val(taxid), path(bundle), path(full_gfa), path(clip_stats), path(identities), path(build_versions)
    path(adapter)
    output:
    tuple val(taxid), path('representation/clipping_report.md'), emit: report
    tuple val(taxid), path('representation/*'), emit: files
    path('representation/versions.tsv'), emit: versions
    script:
    """
    set -euo pipefail
    python3 ${adapter} --bundle ${bundle} --full-gfa ${full_gfa} --clip-stats ${clip_stats} \\
        --identities ${identities} --build-versions ${build_versions} \\
        --image '${task.container}' --out representation
    """
}
