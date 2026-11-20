nextflow.enable.dsl = 2

params.protein_fixture = null
include { CONTAMINATION_PROTEIN_DB } from '../workflows/contamination_protein_db.nf'
include { CONTAMINATION_VERIFY_PROTEIN; CONTAMINATION_PREPARE; CONTAMINATION_BLASTX } from '../modules/contamination_evidence.nf'

workflow {
    def fixture = new groovy.json.JsonSlurper().parseText(file(params.protein_fixture).text)
    def base = file("${projectDir}/..").normalize()
    CONTAMINATION_PROTEIN_DB(fixture.specification + [taxonomy: fixture.taxonomy], base)
    def bundles = CONTAMINATION_PROTEIN_DB.out.directory
    def database = bundles.map { directory -> file(directory).resolve('reference_proteomes.dmnd') }
    def manifest = bundles.map { directory -> file(directory).resolve('protein_database.json') }
    def helpers = ['contamination_review.py', 'contamination_evidence.py'].collect { file(base.resolve('py_scripts/' + it)) }
    CONTAMINATION_VERIFY_PROTEIN(database, manifest, fixture.taxonomy, helpers)
    def inputs = Channel.fromList(['A', 'B']).map { sample -> tuple("fixture_${sample}_hap1", file(fixture.fasta)) }
    CONTAMINATION_PREPARE(inputs, 'post_fcs', 1000, 1000, helpers[0])
    def chunks = CONTAMINATION_PREPARE.out.prepared.flatMap { id, metrics, queries ->
        queries.resolve('chunks.tsv').readLines().drop(1).collect { tuple(id, queries.resolve(it.split('\t')[0])) }
    }
    CONTAMINATION_BLASTX(chunks, database.map { it.toString() }, CONTAMINATION_VERIFY_PROTEIN.out)
    CONTAMINATION_BLASTX.out.toList().view { rows ->
        assert rows.size() == 12: "Missing search chunks: ${rows.size()}"
        assert rows.countBy { it[0].toString() }.values().every { it == 6 }
        def taxids = rows.collectMany { id, hits -> hits.readLines().collect { it.split('\t')[2] } }.toSet()
        assert '7955' in taxids && '83333' in taxids: taxids
        'Managed database and all twelve real DIAMOND searches passed across two assemblies'
    }
}
