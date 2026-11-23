nextflow.enable.dsl = 2

// Scheduling-only fixture: force evaluation of the actual Slurm config closure
// using the same worker input name as the production process. No searches run.
process CONTAMINATION_MEGABLAST_NODE {
    cpus 1
    memory '1 GB'
    input:
    val worker
    output:
    path "worker-${worker}.txt"
    script:
    def hosts = ['crest-c018', 'crest-c021', 'crest-c028', 'crest-c029']
    def placement = task.clusterOptions.toString()
    assert placement == "--exclusive --mem=0 --nodelist=${hosts[(worker as int) - 1]}"
    """
    printf '%s\n' '${placement}' > worker-${worker}.txt
    """
}

process CONTAMINATION_MEGABLAST_COLLECT {
    cpus 1
    memory '1 GB'
    input:
    val worker
    output:
    path "collector-${worker}.txt"
    script:
    def placement = task.clusterOptions.toString()
    def controller = System.getenv('SLURMD_NODENAME') ?: System.getenv('HOSTNAME')
    assert controller && placement == "--nodelist=${controller}"
    """
    printf '%s\n' '${placement}' > collector-${worker}.txt
    """
}

workflow {
    CONTAMINATION_MEGABLAST_NODE(Channel.of(1, 2, 3, 4))
    CONTAMINATION_MEGABLAST_COLLECT(Channel.of(1, 2, 3, 4))
    CONTAMINATION_MEGABLAST_NODE.out.toList().view { rows ->
        assert rows.size() == 4
        'All four exclusive worker placements resolved correctly'
    }
}
