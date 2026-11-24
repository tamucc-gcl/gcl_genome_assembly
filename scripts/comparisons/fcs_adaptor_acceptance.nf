nextflow.enable.dsl=2
include { FCS_ADAPTOR as SCREEN; FCS_ADAPTOR_CLEAN as CLEAN; FCS_ADAPTOR as RESCREEN } from '../../modules/fcs_adaptor.nf'

process FIXTURE {
    conda 'conda-forge::python=3.11'
    cpus 1
    memory '1 GB'
    input:
    path helpers
    output:
    path 'fixture.fasta'
    script:
    """
    python3 fcs_adaptor_acceptance.py fetch --out fixture.fasta
    """
}
process CHECK {
    conda 'conda-forge::python=3.11'
    cpus 1
    memory '1 GB'
    publishDir params.outdir,mode:'copy'
    input:
    path source, stageAs: 'source/*'
    path report, stageAs: 'initial/*'
    path clean, stageAs: 'cleaned/*'
    path rescreen, stageAs: 'rescreen/*'
    path helpers
    output:
    path 'acceptance.json'
    script:
    """
    python3 fcs_adaptor_acceptance.py check --source ${source} --report ${report} \
        --clean ${clean} --rescreen ${rescreen} --out acceptance.json
    """
}
workflow {
    def base=file("${projectDir}/../..")
    def helpers=['tests/fcs_adaptor_acceptance.py','py_scripts/contamination_review.py',
        'py_scripts/contamination_cohort_report.py'].collect { file(base.resolve(it),checkIfExists:true) }
    FIXTURE(helpers)
    def launcher=file(base.resolve('py_scripts/run_fcs_adaptor.py'),checkIfExists:true)
    SCREEN(FIXTURE.out.map { fa -> tuple([id:'ncbi_fixture'],fa,'euk','nextflow','contig') },launcher)
    CLEAN(FIXTURE.out.combine(SCREEN.out.action_report).map { fa,m,report -> tuple(m,fa,report,'contig') })
    RESCREEN(CLEAN.out.cleaned_fasta.map { m,fa -> tuple(m,fa,'euk','nextflow','final') },launcher)
    CHECK(FIXTURE.out,SCREEN.out.action_report.map { m,p -> p },CLEAN.out.cleaned_fasta.map { m,p -> p },
        RESCREEN.out.action_report.map { m,p -> p },helpers)
}
