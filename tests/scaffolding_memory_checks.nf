nextflow.enable.dsl=2
include { SCAFFOLD_HIC } from '__REPO__/modules/scaffold_hic.nf'
workflow {
    def root=file('/tmp/chimera-nxf-runtime/scaffolding-fixtures');root.mkdirs()
    def fasta=file("${root}/input.fa");fasta.text='>s\nACGT\n'
    def links=file("${root}/input.bin");links.text='fixture'
    def index=file("${root}/unused_index");index.text='unused for BIN input'
    SCAFFOLD_HIC(Channel.of([[id:'fixture'],links,index,fasta,'round1',
        [min_contig_length:20000,resolutions:'20000,50000',min_mapq:10]]))
    SCAFFOLD_HIC.out.agp.view { meta,agp ->
        assert agp.text.contains('s\t1\t4')
        'PASS: production scaffolding completed with requested-resolution memory policy'
    }
}
