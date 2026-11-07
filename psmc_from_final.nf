// Removable comparison entry: the production PSMC workflow on accepted final inputs.
nextflow.enable.dsl=2
include { PSMC_ANALYSIS } from './workflows/psmc.nf'
workflow {
    if(params.psmc_mode!='assemblies') error 'This bounded comparison entry requires psmc_mode=assemblies'
    PSMC_ANALYSIS(Channel.value(file(params.psmc_final_eligibility,checkIfExists:true)),Channel.empty())
}
