/* Check the pinned image before allocating graph-build resources. */
process PANGENOME_PREFLIGHT {
    label 'cactus_preflight'
    input:
    val requested
    output:
    path 'cactus_capabilities.txt', emit: ready
    script:
    """
    set -euo pipefail
    cactus-pangenome --help > cactus_capabilities.txt 2>&1
    for flag in --gref --haplo --refContigs --vcfbub; do
        if ! grep -q -- "\$flag" cactus_capabilities.txt; then
            echo "Configured Cactus image lacks required option \$flag; no graph was started." >&2
            exit 1
        fi
    done
    """
}