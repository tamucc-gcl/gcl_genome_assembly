/*
========================================================================================
    SETUP_CHLOE MODULE
========================================================================================
    Repo location: modules/setup_chloe.nf

    Installs Chloë (angiosperm plastid annotator) once into a permanent shared dir, before
    any parallel ANNOTATE_PLASTID runs. Chloë is a Julia tool: it needs its code repo and a
    references repo cloned as siblings, plus the Julia project instantiated. We treat that
    like a database — set up once, guard with a sentinel, reuse everywhere.

    Layout created under <chloe_dir>:
      chloe/              (github.com/ian-small/chloe — the code, incl. chloe.jl)
      chloe_references/   (github.com/ian-small/chloe_references — auto-found as a sibling)
      depot/              (JULIA_DEPOT_PATH — the instantiated packages, shared across nodes)

    Input:
      - chloe_dir:  absolute dir (params.chloe_downloads)
      - force:      re-setup even if the sentinel exists
    Output:
      - dir:        chloe_dir (val), consumed by ANNOTATE_PLASTID
========================================================================================
*/

process SETUP_CHLOE {
    tag "chloe_setup"
    label 'chloe'

    input:
    val chloe_dir
    val force

    output:
    val "${chloe_dir}", emit: dir

    script:
    """
    set -eu
    mkdir -p "${chloe_dir}"
    export JULIA_DEPOT_PATH="${chloe_dir}/depot"
    SENTINEL="${chloe_dir}/.chloe_ready"

    if [ "${force}" = "true" ] || [ ! -f "\${SENTINEL}" ]; then
        [ -d "${chloe_dir}/chloe" ]            || git clone --depth 1 https://github.com/ian-small/chloe            "${chloe_dir}/chloe"
        [ -d "${chloe_dir}/chloe_references" ] || git clone --depth 1 https://github.com/ian-small/chloe_references "${chloe_dir}/chloe_references"
        ( cd "${chloe_dir}/chloe" && julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' )
        date -Is > "\${SENTINEL}"
        echo "[SETUP_CHLOE] Chloë ready at ${chloe_dir}"
    else
        echo "[SETUP_CHLOE] present (sentinel found) — skipping"
    fi
    """

    stub:
    """
    mkdir -p "${chloe_dir}/chloe" "${chloe_dir}/chloe_references" "${chloe_dir}/depot"
    touch "${chloe_dir}/.chloe_ready"
    """
}
