/* Cheap current-run status: an absent exact join test is not a clean assembly result. */
process CHIMERA_COORDINATE_SUMMARY {
    cpus 1
    memory '1 GB'
    time '10m'
    publishDir "${params.outdir}/assembly/chimeras", mode: params.publish_dir_mode
    input:
    val assemblies
    val tested_ids
    val older_status
    output:
    path 'chimera_coordinate_status.md', emit: report
    script:
    // Normalize GString metadata IDs and process-output String IDs before hash lookup.
    def tested = tested_ids.collect { it.toString() } as Set
    def clean = { x -> x.toString().replaceAll(/[\r\n\t|]/, ' ') }
    def older = older_status.collectEntries { r -> [(r.assembly.toString()): r] }
    def recoveryRows = assemblies.sort { it.id }.collect { r ->
        def item = older[r.id.toString()]
        def recovered = item ? (item.statuses.recovered_exact_flanks_and_gap ?: 0) : null
        def status = item ? "${recovered}/${item.total} recovered; ${item.total - recovered} unresolved; ${item.transition_intervals ?: 0} transition intervals; ${item.diagnostic_profiles ?: 0} diagnostic profiles requested" : 'Not run'
        "| ${clean(r.id)} | ${status} |"
    }.join('\n')
    def rows = assemblies.sort { it.id }.collect { r ->
        def state = tested.contains(r.id.toString()) ? 'AGP/FASTA checked; last-round join scope' :
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
Older joins are assessed separately using exact mapped flanks and a retained N-gap.
Recovery is conservative: altered, missing, ambiguous or split flanks stay unresolved.
A recovered coordinate is not proof of a biological misjoin or permission to cut.
No older join is automatically cut by this recovery branch.

| Assembly | Older-join recovery |
|---|---|
${recoveryRows}

Detailed audits, review tables and *.transition_intervals.tsv are under assembly/chimeras/older_joins.
Interval midpoints are diagnostic positions, not safe cuts. Each evidence TSV records figure_status.
Unresolved older joins are **not cleared**.
Current-assembly harmonization can still flag those scaffolds for review.
These evidence coordinates describe the pre-finishing assembly, not the final FASTA.
A checked coordinate frame does not mean that Hi-C supports a cut or that the assembly is chimera-free.

</details>
COORDINATES
    """
}