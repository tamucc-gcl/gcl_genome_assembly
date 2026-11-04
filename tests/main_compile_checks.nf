nextflow.enable.dsl=2
include { VALIDATE_INPUTS } from '__REPO__/main.nf'
workflow {
    VALIDATE_INPUTS()
}
