nextflow.enable.dsl = 2

include { assemblyKeys; hitsByAssembly } from '../workflows/contamination_search.nf'

process CHECK_ASSEMBLY {
    tag "${id}"
    input:
    tuple val(id), val(proteinCount), val(nucleotideCount), val(expected)
    output:
    tuple val(id), path('verified.txt')
    script:
    assert proteinCount == expected
    assert nucleotideCount == expected
    """
    printf '%s\\n' '${id}: ${expected} protein and nucleotide chunks' > verified.txt
    """
}

workflow {
    // The production cohort has 2,088 chunks across ten haplotypes.
    def sizes = [218, 218, 200, 209, 215, 213, 212, 194, 193, 216, 0]
    def specifications = sizes.withIndex().collect { size, index ->
        [id: "sample_${index}_hap1", size: size]
    }
    specifications.each { spec ->
        assert nextflow.util.CacheHelper.hasher(spec.id).hash() ==
            nextflow.util.CacheHelper.hasher(spec.id.toString()).hash()
    }
    assert sizes.sum() == 2088
    def counts = Channel.fromList(specifications.collect { tuple(it.id, it.size) })
    // Process/cache outputs can carry plain strings while upstream metadata uses GStrings.
    def proteinHits = Channel.fromList(specifications.collectMany { spec ->
        (0..<spec.size).collect { tuple(spec.id.toString(), file("${spec.id}/protein_${it}.tsv")) }
    })
    def nucleotideHits = Channel.fromList(specifications.reverse().collectMany { spec ->
        (0..<spec.size).collect { tuple(spec.id, file("${spec.id}/nucleotide_${it}.tsv")) }
    })
    def sources = assemblyKeys(Channel.fromList(specifications.collect { tuple(it.id, it.size) }))
    def prepared = assemblyKeys(Channel.fromList(specifications.collect { tuple(it.id.toString(), 'prepared') }))
    def coverage = assemblyKeys(Channel.fromList(specifications.collect { tuple(it.id, [], [], []) }))
    def input = sources.join(prepared, failOnMismatch: true, failOnDuplicate: true)
        .join(hitsByAssembly(counts, proteinHits), failOnMismatch: true, failOnDuplicate: true)
        .join(hitsByAssembly(counts, nucleotideHits), failOnMismatch: true, failOnDuplicate: true)
        .join(coverage, failOnMismatch: true, failOnDuplicate: true)
        .map { id, expected, preparedDir, proteins, nucleotides, bam, bai, cov ->
            assert id instanceof String
            assert bam == [] && bai == [] && cov == []
            assert proteins*.name.unique().size() == expected
            assert nucleotides*.name.unique().size() == expected
            tuple(id, proteins.size(), nucleotides.size(), expected)
        }
    CHECK_ASSEMBLY(input)
    CHECK_ASSEMBLY.out.toList().view { rows ->
        assert rows.size() == specifications.size()
        'Assembly grouping passed: 2,088 chunks per evidence source, ten haplotypes plus an empty assembly'
    }
}
