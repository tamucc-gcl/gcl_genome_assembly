// Removable diagnostic: run the production module once using retained Hi-C BIN.
nextflow.enable.dsl=2
include { SCAFFOLD_HIC } from '../../modules/scaffold_hic.nf'
workflow {
    def fasta=file(params.replay_fasta,checkIfExists:true)
    def links=file(params.replay_links,checkIfExists:true)
    def expected=file(params.replay_expected_agp,checkIfExists:true)
    def unusedIndex=file("${projectDir}/../../assets/NO_PAIRS",checkIfExists:true)
    def meta=[id:params.replay_assembly]
    def settings=[min_contig_length:params.yahs_round1_min_contig_bp,
        min_mapq:params.yahs_round1_min_mapq,resolutions:params.yahs_round1_resolutions,
        rounds_per_resolution:params.yahs_round1_rounds_per_resolution,
        enzyme:params.yahs_round1_enzyme,
        no_contig_ec:!params.yahs_round1_contig_ec,
        no_scaffold_ec:!params.yahs_round1_scaffold_ec]
    SCAFFOLD_HIC(Channel.of([meta,links,unusedIndex,fasta,'round1',settings]))
    SCAFFOLD_HIC.out.agp.view { identity,agp ->
        def matching=agp.text==expected.text
        def result=[assembly:identity.id,exact_agp_match:matching,
                    expected_agp:expected.toString(),replayed_agp:agp.toString()]
        file("${params.outdir}/../replay_check.json").text=
            groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(result))+'\n'
        assert matching : 'Replayed scaffolding differs: inspect replay_check.json and YaHS log'
        'PASS: production YaHS replay exactly matches accepted AGP'
    }
}
