/* Cheap current-run status: an absent exact join test is not a clean assembly result. */
process CHIMERA_COORDINATE_SUMMARY {
    cpus 1
    memory '1 GB'
    time '10m'
    publishDir "${params.outdir}/assembly/chimeras", mode: params.publish_dir_mode
    input:
    val assemblies
    val tested_ids
    output:
    path 'chimera_coordinate_status.md', emit: report
    script:
    def tested = tested_ids as Set
    def clean = { x -> x.toString().replaceAll(/[\r\n\t|]/, ' ') }
    def rows = assemblies.sort { it.id }.collect { r ->
        def state = tested.contains(r.id) ? 'AGP/FASTA checked; last-round join scope' :
            (!params.chimera_detect ? 'Detection disabled' :
             (!r.hic ? 'No Hi-C scaffold-join evidence' :
              (!r.harmonized ? 'No harmonization reference/map' : 'Exact join evidence unavailable')))
        "| ${clean(r.id)} | ${state} |"
    }.join('\n')
    """
    cat > chimera_coordinate_status.md <<'COORDINATES'
<details>
<summary>Chimera assessment coordinates and coverage</summary>

Assessment and optional cuts precede gap filling, telomere extension, and final renaming.
Hi-C pairs are projected through the AGP of the last scaffolding round from mappings
to that round's input. No projection is made through Inspector corrections.

| Assembly | Coordinate/evidence status |
|---|---|
${rows}

Only joins recorded by the last scaffolding round are tested for exact cuts.
Older joins inside corrected components have unresolved provenance and are **not cleared**.
Current-assembly harmonization can still flag those scaffolds for review.
These evidence coordinates describe the pre-finishing assembly, not the final FASTA.
A checked coordinate frame does not mean that Hi-C supports a cut or that the assembly is chimera-free.

</details>
COORDINATES
    """
}