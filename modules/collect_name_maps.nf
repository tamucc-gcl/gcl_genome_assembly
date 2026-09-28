/* Collect final per-assembly name maps into assembly/old_name/new_name columns for report labels. */

process COLLECT_NAME_MAPS {
    tag "name_maps"
    label 'process_single'

    publishDir "${params.outdir}/assembly/harmonization", mode: params.publish_dir_mode

    input:
    path(name_maps)

    output:
    path("all_name_maps.tsv"), emit: map

    script:
    """
    set -euo pipefail
    printf 'assembly\\told_name\\tnew_name\\n' > all_name_maps.tsv
    for f in ${name_maps}; do
        id=\$(basename "\$f" .name_map.tsv)
        # per-assembly map header: old_name  new_name  orient  order  length  class  ref_span  flags
        tail -n +2 "\$f" | awk -v a="\$id" 'BEGIN{FS=OFS="\\t"} NF>=2 {print a, \$1, \$2}' >> all_name_maps.tsv
    done
    """

    stub:
    """
    printf 'assembly\\told_name\\tnew_name\\n' > all_name_maps.tsv
    """
}
