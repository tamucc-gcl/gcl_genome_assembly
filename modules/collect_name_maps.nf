/*
========================================================================================
    COLLECT NAME MAPS MODULE
========================================================================================
    Repo location: modules/collect_name_maps.nf

    Merges the per-assembly harmonization name-maps (<id>.harmonized_name_map.tsv, emitted
    by HARMONIZE_SCAFFOLDS as the 3rd element of its `assemblies` tuple) into a single
    3-column table keyed by assembly id, for the teloclip contig remap in the summary report
    (pre-harmonization scaffold_N -> final chrN_1 / unplaced_N).

    Input : path(name_maps)   one or more <id>.harmonized_name_map.tsv
    Output: map = all_name_maps.tsv   (columns: assembly  old_name  new_name)
========================================================================================
*/

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
        id=\$(basename "\$f" .harmonized_name_map.tsv)
        # per-assembly map header: old_name  new_name  orient  order  length  class  ref_span  flags
        tail -n +2 "\$f" | awk -v a="\$id" 'BEGIN{FS=OFS="\\t"} NF>=2 {print a, \$1, \$2}' >> all_name_maps.tsv
    done
    """

    stub:
    """
    printf 'assembly\\told_name\\tnew_name\\n' > all_name_maps.tsv
    """
}
