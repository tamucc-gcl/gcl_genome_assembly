nextflow.enable.dsl=2
include { PSMC_ANALYSIS } from '__REPO__/workflows/psmc.nf'
workflow {
    def root=file(params.psmc_fixture)
    PSMC_ANALYSIS(Channel.value(file("${root}/eligibility.json")),
        Channel.value(tuple('1',file("${root}/reference.json"))),
        Channel.value(tuple([sample:'read_ind'],file("${root}/complete.fastq"))))
    PSMC_ANALYSIS.out.report.view { p ->
        if(p.name=='summary.tsv') {
            if(params.psmc_test_case!='all_skipped') {
                def result=params.psmc_test_case=='poor_callability' ? 'not_fitted' : 'completed'
                assert p.text.contains("assembly_ind\t1\tassemblies\t${result}")
                assert p.text.contains("read_ind\t1\thifi\t${result}")
            }
            assert p.text.contains('haploid_ind\t1\tnot_run\tskipped')
        }
        "PSMC product: ${p}"
    }
}
