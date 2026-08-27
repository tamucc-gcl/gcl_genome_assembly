/*
========================================================================================
    PANGENOME REARRANGE MODULE
========================================================================================
    Repo location: modules/pangenome_rearrange.nf

    Collects every per-chromosome PANGENOME_UNTANGLE table for one graph flavour and turns
    them into rearrangement calls plus a CANDIDATE TABLE.

    The candidate table is the reason this module exists rather than the BED being published
    directly. A segregating inversion polymorphism on chr10 was found during development --
    ~12 Mb of reference footprint, carried by 2 of 10 haplotypes, both carriers' sister
    haplotypes clean, so both individuals heterozygous at AF 0.2 -- and finding it required
    several rounds of ad-hoc awk over raw untangle output. It appeared in no published
    artifact. This module makes it the top row of a table, and that is the acceptance test
    for the reporting layer: if a result like that is not visible without hand analysis, the
    reporting is incomplete.

    The artifact discriminator is a join against harmonization's own flags rather than any
    projection statistic. A scaffold spanning two reference chromosomes has no correct single
    orientation, so one half projects as 100% inverted -- but harmonization already knows,
    because `unsupported(ref4+ref10:0f/8s)` records that the junction was seen fused in 0
    voters and split in 8. Zero fused voters beats "100% inverted" as evidence. Both flag
    spellings are matched: `unsupported(...)` and `chimera_suspect(...)`.

    The harmonization report is OPTIONAL. Without it the tables are still produced and the
    audit records `harmonization_rows_loaded=0`, so a missing join is visible rather than
    silently degrading the artifact flags.

    Input : tuple(taxid, flavor, tsvs, harmonization_report)  -- the report is JOINED INTO
            the tuple rather than passed as a separate channel. Passing it separately emits
            one item while the grouped untangle channel emits one per flavour, so Nextflow
            would run this process once and SILENTLY DROP the second flavour.
    Output: candidates / orientation / inversions BED / duplications / audit / versions
========================================================================================
*/

process PANGENOME_REARRANGE {
    tag "${taxid}:${flavor}"
    label 'pangenome_rearrange'

    publishDir "${params.outdir}/pangenome/${taxid}/rearrange", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), path(tsvs), path(harmonization)
    path(script)

    output:
    tuple val(taxid), val(flavor),
          path("${taxid}.${flavor}.rearrangement_candidates.tsv"), emit: candidates
    tuple val(taxid), val(flavor),
          path("${taxid}.${flavor}.query_orientation.tsv"),        emit: orientation
    tuple val(taxid), val(flavor),
          path("${taxid}.${flavor}.inversions.bed"),               emit: inversions
    tuple val(taxid), val(flavor),
          path("${taxid}.${flavor}.duplications.tsv"),             emit: duplications
    tuple val(taxid), val(flavor),
          path("${taxid}.${flavor}.untangle_audit.tsv"),           emit: audit
    path("versions.tsv"),                                          emit: versions

    script:
    def harm = (harmonization && harmonization.name != 'NO_FILE')
                 ? "--harmonization ${harmonization}" : ''
    def minjac = params.pangenome_untangle_min_jaccard ?: 0.1
    def minseg = params.pangenome_inv_min_bp           ?: 1000
    def gap    = params.pangenome_inv_run_merge_gap    ?: 200000
    def afrac  = params.pangenome_untangle_orientation_artifact_frac ?: 0.95
    """
    set -euo pipefail

    python3 ${script} \\
        --untangle ${tsvs} \\
        --outdir . \\
        --label ${taxid} \\
        --flavor ${flavor} \\
        --min-jaccard ${minjac} \\
        --min-seg-bp ${minseg} \\
        --run-merge-gap ${gap} \\
        --artifact-frac ${afrac} \\
        ${harm}

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    for f in rearrangement_candidates query_orientation duplications untangle_audit; do
        : > ${taxid}.${flavor}.\$f.tsv
    done
    : > ${taxid}.${flavor}.inversions.bed
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
