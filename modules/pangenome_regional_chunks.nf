process PANGENOME_REGIONAL_CHUNK_ASSIGN {
    tag "${taxid}:${policy}:${variant}:${chunk}"
    label 'cactus_tools'
    cpus 1
    memory '16 GB'
    time '2h'
    publishDir "${params.outdir}/pangenome/${taxid}/deletion_filter", mode: params.publish_dir_mode,
        pattern: 'chunk_*'
    input:
    tuple val(taxid), val(policy), path(gfa), path(paf), path(config), path(selection), path(identities),
        path(provenance), val(contigs), path(binary), path(build_info), path(tools), val(variant), val(chunk)
    path(adapter)
    path(helpers)
    output:
    tuple val(taxid), val(policy), val(variant), val(chunk), path("chunk_${policy}_${variant}_${chunk}"), emit: results
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
    python3 ${adapter} run --prepared prepared --tools ${tools} --fixed-binary ${binary} --fixed-info ${build_info} \\
        --variant ${variant} --chunk ${chunk} --contigs ${contigs.join(' ')} --out chunk_${policy}_${variant}_${chunk}
    """
}

process PANGENOME_REGIONAL_CHUNK_REPORT {
    tag "${taxid}"
    label 'cactus_tools'
    cpus 1
    memory '4 GB'
    time '30m'
    publishDir "${params.outdir}/pangenome/${taxid}/deletion_filter", mode: params.publish_dir_mode,
        pattern: 'chunk_summary/*', saveAs: { it.replaceFirst('^chunk_summary/', '') }
    input:
    tuple val(taxid), path(baselines), path(comparisons)
    path(adapter)
    path(helpers)
    output:
    tuple val(taxid), path('chunk_summary/regional_chunk_report.md'), emit: report
    tuple val(taxid), path('chunk_summary/*'), emit: files
    script:
    """
    set -euo pipefail
    python3 ${adapter} report --baselines ${baselines} --comparisons ${comparisons} --out chunk_summary
    """
}

workflow PANGENOME_REGIONAL_CHUNK_TEST {
    take:
    filtered
    fixed_tool
    trace_tools
    baselines
    adapter
    helpers
    main:
    ch_cases = filtered.combine(fixed_tool).combine(trace_tools).flatMap { row ->
        ['native','fixed'].collectMany { variant -> [100000,25000].collect { chunk -> row + [variant,chunk] } }
    }
    PANGENOME_REGIONAL_CHUNK_ASSIGN(ch_cases, adapter, helpers)
    PANGENOME_REGIONAL_CHUNK_REPORT(baselines.map { taxid, policy, dir -> tuple(taxid,dir) }.groupTuple()
        .join(PANGENOME_REGIONAL_CHUNK_ASSIGN.out.results.map { taxid, policy, variant, chunk, dir -> tuple(taxid,dir) }.groupTuple()),
        adapter, helpers)
    emit:
    report = PANGENOME_REGIONAL_CHUNK_REPORT.out.report
}
