/* Optional diagnostic batch. All construction products remain read-only. */
process PANGENOME_ASSIGNMENT_PREPARE {
    tag "${taxid}"
    label 'cactus_tools'
    cpus 1
    memory '8 GB'
    time '2h'
    publishDir "${params.outdir}/pangenome/${taxid}/assignment/preparation", mode: params.publish_dir_mode,
        pattern: 'prepared/*', saveAs: { fn ->
            fn.endsWith('.json') || fn.endsWith('.xml') || fn.endsWith('.tsv') ? fn.replaceFirst('^prepared/', '') : null
        }
    input:
    tuple val(taxid), path(gfa), path(paf), val(contigs), path(identities), path(sources)
    path(adapters)
    output:
    tuple val(taxid), path('prepared/cactus.sv.gfa'), path('prepared/filtered.paf'),
        path('prepared/effective_config.xml'), path('prepared/selection.tsv'),
        path('prepared/identities.tsv'), path('prepared/preparation.json'), val(contigs), emit: prepared
    script:
    """
    set -euo pipefail
    python3 pangenome_assignment_replay.py prepare --gfa ${gfa} --paf ${paf} \\
        --identities ${identities} --sources ${sources} --out prepared
    """
}

process PANGENOME_ASSIGNMENT_REPLAY {
    tag "${taxid}:${case_name}"
    label 'cactus_tools'
    cpus 1
    memory '16 GB'
    time '2h'
    publishDir "${params.outdir}/pangenome/${taxid}/assignment", mode: params.publish_dir_mode,
        pattern: 'replay_*'
    input:
    tuple val(taxid), path(gfa), path(paf), path(config), path(selection),
        path(identities), path(provenance), val(contigs), val(case_name)
    path(adapters)
    output:
    tuple val(taxid), val(case_name), path("replay_${case_name}"), emit: results
    script:
    """
    set -euo pipefail
    mkdir prepared
    ln -s "\$PWD/${gfa}" prepared/cactus.sv.gfa
    ln -s "\$PWD/${paf}" prepared/filtered.paf
    ln -s "\$PWD/${config}" prepared/effective_config.xml
    ln -s "\$PWD/${selection}" prepared/selection.tsv
    ln -s "\$PWD/${identities}" prepared/identities.tsv
    ln -s "\$PWD/${provenance}" prepared/preparation.json
    python3 pangenome_assignment_replay.py replay --prepared prepared --case '${case_name}' \\
        --contigs ${contigs.join(' ')} --out replay_${case_name}
    """
}

process PANGENOME_ASSIGNMENT_REPORT {
    tag "${taxid}"
    label 'cactus_tools'
    cpus 1
    memory '1 GB'
    time '10m'
    publishDir "${params.outdir}/pangenome/${taxid}/assignment", mode: params.publish_dir_mode,
        pattern: 'summary/*', saveAs: { it.replaceFirst('^summary/', '') }
    input:
    tuple val(taxid), path(directories)
    path(adapters)
    output:
    tuple val(taxid), path('summary/assignment_report.md'), emit: report
    tuple val(taxid), path('summary/*'), emit: files
    script:
    """
    set -euo pipefail
    python3 pangenome_assignment_replay.py report --directories ${directories} --out summary
    """
}
