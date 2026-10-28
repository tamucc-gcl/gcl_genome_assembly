nextflow.enable.dsl=2
include { INPUT_TAXONOMY } from './workflows/input_taxonomy.nf'
include { DOWNLOAD_TAXDUMP } from './modules/download_taxdump.nf'
include { RESOLVE_TAXONOMY } from './modules/resolve_taxonomy.nf'
include { FIND_MITO_REFERENCE } from './modules/find_mito_reference.nf'

process IDENTITY_REFERENCE_EVIDENCE {
    cpus 1
    memory '1 GB'
    time '10m'
    publishDir "${params.outdir}", mode: params.publish_dir_mode
    input:
    tuple val(taxid), path(taxonomy)
    tuple val(reference_taxid), path(reference_info)
    output:
    path 'resolved_taxonomy.tsv'
    path 'mitochondrial_reference.tsv'
    path 'identity_reference_acceptance.md'
    script:
    assert taxid.toString()==reference_taxid.toString()
    """
    set -euo pipefail
    cp '${taxonomy}' resolved_taxonomy.tsv
    cp '${reference_info}' mitochondrial_reference.tsv
    cat > identity_reference_acceptance.md <<'ACCEPTANCE'
# Native identity/reference check

PASS: exact scientific-name resolution agrees with the supplied taxid. Both
sample identities survived, the installed taxonkit resolved the lineage, and
findMitoReference produced a nonempty, matching reference pair passing the
configured minimum length and organism checks.

This tests installed tools and reference access, not an assembly or graph rebuild.
The reference may be from a related species. It is not independently certified
as the best available biological reference by this check.

See resolved_taxonomy.tsv and mitochondrial_reference.tsv for the actual records.
ACCEPTANCE
    """
}

workflow {
    if (!params.taxid || !params.species) error 'Supply --taxid and --species for the native equivalence check'
    DOWNLOAD_TAXDUMP(Channel.value(file(params.diamond_taxdump_dir)),Channel.value(false))
    INPUT_TAXONOMY(Channel.of(
        tuple([sample:'by_taxid',taxid:params.taxid.toString(),species:null],[:]),
        tuple([sample:'by_name',taxid:null,species:params.species.toString()],[:])),
        DOWNLOAD_TAXDUMP.out.taxdump_dir)
    resolved=INPUT_TAXONOMY.out.samples.toList().map { rows ->
        assert rows*.getAt(0)*.sample.sort()==['by_name','by_taxid']
        assert rows.every { it[0].taxid.toString()==params.taxid.toString() }
        params.taxid.toString()
    }
    RESOLVE_TAXONOMY(resolved.combine(DOWNLOAD_TAXDUMP.out.taxdump_dir))
    FIND_MITO_REFERENCE(resolved.map { taxid -> tuple(taxid,params.species.toString()) })
    IDENTITY_REFERENCE_EVIDENCE(RESOLVE_TAXONOMY.out.tsv,FIND_MITO_REFERENCE.out.ref_info)
}
