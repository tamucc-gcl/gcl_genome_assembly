/*
========================================================================================
    RESOLVE_TAXONOMY — taxid -> name + ranks + full lineage (local taxdump, offline)
========================================================================================
    Once per DISTINCT taxid. main.nf parses the TSV and applies functions/taxonomy.nf.
========================================================================================
*/
process RESOLVE_TAXONOMY {
    tag "${taxid}"
    label 'taxonkit'

    input:
    tuple val(taxid), path(taxdump_dir)

    output:
    tuple val(taxid), path("${taxid}.taxonomy.tsv"), emit: tsv

    script:
    """
    set -euo pipefail
    if [ ! -s "${taxdump_dir}/names.dmp" ] || [ ! -s "${taxdump_dir}/nodes.dmp" ]; then
        echo "[RESOLVE_TAXONOMY] ERROR: names.dmp/nodes.dmp not in ${taxdump_dir}" >&2; exit 1
    fi

    # Ranks. Codes: {k}=superkingdom {K}=kingdom {p}=phylum {c}=class {o}=order {f}=family {g}=genus {s}=species
    RANKS=\$(echo "${taxid}" \\
        | taxonkit reformat -I 1 --data-dir "${taxdump_dir}" -F -r NA \\
              -f "{k}\\t{K}\\t{p}\\t{c}\\t{o}\\t{f}\\t{g}\\t{s}" \\
        | cut -f2-)

    # Full lineage (semicolon-joined clade NAMES) — captures unranked clades NCBI hides from
    # {c}: Liliopsida (monocots), eudicotyledons, Endopterygota, etc. The BUSCO map matches this.
    LINEAGE=\$(echo "${taxid}" | taxonkit lineage --data-dir "${taxdump_dir}" | cut -f2)

    printf 'taxid\\tsuperkingdom\\tkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\tlineage\\n' > ${taxid}.taxonomy.tsv
    printf '%s\\t%s\\t%s\\n' "${taxid}" "\${RANKS}" "\${LINEAGE}" >> ${taxid}.taxonomy.tsv
    """

    stub:
    """
    printf 'taxid\\tsuperkingdom\\tkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\tlineage\\n' > ${taxid}.taxonomy.tsv
    printf '${taxid}\\tEukaryota\\tMetazoa\\tChordata\\tActinopteri\\tClupeiformes\\tClupeidae\\tSpratelloides\\tSpratelloides delicatulus\\tcellular organisms;Eukaryota;Metazoa;Chordata;Actinopteri;Clupeiformes;Clupeidae;Spratelloides;Spratelloides delicatulus\\n' >> ${taxid}.taxonomy.tsv
    """
}