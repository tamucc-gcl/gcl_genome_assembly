nextflow.enable.dsl=2
params.run_pangenome = (System.getenv('CORE_INPUT_CASE') == 'graph_with_invalid_sample')
params.hifiasm_primary = false
params.outdir = '/tmp/chimera-nxf-runtime/reference-results'
params.publish_dir_mode = 'copy'
include { parseSampleSheet } from '__REPO__/functions/input_validation.nf'
include { SHARED_REFERENCE } from '__REPO__/modules/shared_reference.nf'

workflow {
    def root = file('/tmp/chimera-nxf-runtime/input-checks')
    root.mkdirs()
    ['hifi.bam', 'one.fq', 'two.fq', 'hic1.fq', 'hic2.fq', 'hic3.fq', 'hic4.fq'].each {
        file("${root}/${it}").text = 'fixture\n'
    }
    def sheet = file("${root}/samples.csv")
    def hic = file("${root}/hic.csv")
    def check = { text, hicText, expected, invalid, skipped ->
        sheet.text = text
        hic.text = hicText
        def result = parseSampleSheet(sheet.toString(), hicText ? hic.toString() : null)
        assert result.samples*.getAt(0)*.sample == expected
        assert result.status.findAll { it.status == 'invalid' }.size() == invalid
        assert result.status.findAll { it.status == 'skipped' }.size() == skipped
        result
    }
    def header = 'sample_id,hifi_bam,sr_r1,sr_r2,hic_r1,hic_r2,ploidy,mito_tool\n'
    def negative = System.getenv('CORE_INPUT_CASE') ?: 'success'
    if (negative == 'duplicate_sample')
        check(header + 'same,hifi.bam,,,,,diploid,none\nsame,,,,,,,\n', '', [], 0, 0)
    if (negative == 'duplicate_file')
        check(header + 'a,hifi.bam,,,,,diploid,none\nb,hifi.bam,,,,,diploid,none\n', '', [], 0, 0)
    if (negative == 'unknown_hic_sample')
        check(header + 'a,hifi.bam,,,,,diploid,none\n',
            'sample_id,library_id,readset_id,hic_r1,hic_r2\nabsent,A,A1,hic1.fq,hic2.fq\n', [], 0, 0)
    if (negative == 'no_valid_samples')
        check(header + 'a,,,,hic1.fq,hic2.fq,diploid,none\n', '', [], 0, 1)
    if (negative == 'duplicate_header')
        check('sample_id,hifi_bam,hifi_bam\na,hifi.bam,hifi.bam\n', '', [], 0, 0)
    if (negative == 'malformed_record')
        check('sample_id,hifi_bam\na,hifi.bam,unexpected\n', '', [], 0, 0)
    if (negative == 'unsafe_sample')
        check(header + '../a,hifi.bam,,,,,diploid,none\n', '', [], 0, 0)
    if (negative == 'duplicate_readset')
        check(header + 'a,hifi.bam,,,,,diploid,none\n',
            'sample_id,library_id,readset_id,hic_r1,hic_r2\na,A,R,hic1.fq,hic2.fq\na,B,R,hic3.fq,hic4.fq\n', [], 0, 0)
    if (negative == 'missing_hic_mate')
        check(header + 'a,hifi.bam,,,,,diploid,none\n',
            'sample_id,library_id,readset_id,hic_r1,hic_r2\na,A,R,hic1.fq,\n', [], 0, 0)
    if (negative == 'conflicting_hic_sources')
        check(header + 'a,hifi.bam,,,hic1.fq,hic2.fq,diploid,none\n',
            'sample_id,library_id,readset_id,hic_r1,hic_r2\na,A,R,hic3.fq,hic4.fq\n', [], 0, 0)
    if (negative == 'graph_with_invalid_sample') {
        check(header + 'good,hifi.bam,,,,,diploid,none\nbad,,missing.fq,two.fq,,,diploid,none\n', '', [], 1, 0)
    }
    def mixed = check(header + 'long,hifi.bam,,,,,diploid,none\nshort,,one.fq,two.fq,,,diploid,none\nhic_only,,,,hic1.fq,hic2.fq,diploid,none\n',
        '', ['long','short'], 0, 1)
    assert mixed.samples[0][0].n_hap == 2
    assert mixed.samples[1][0].n_hap == 1
    assert mixed.capabilities.hifiasm && mixed.capabilities.spades
    def libraries = 'sample_id,library_id,readset_id,hic_r1,hic_r2\nlong,A,A1,hic1.fq,hic2.fq\nlong,B,B1,hic3.fq,hic4.fq\n'
    def two = check(header + 'long,hifi.bam,,,,,diploid,none\n', libraries, ['long'], 0, 0)
    assert two.samples[0][1].hic_sets.size() == 2
    assert two.samples[0][1].hic_sets*.library_id == ['A','B']
    def shortHic = check(header + 'short,,one.fq,two.fq,hic1.fq,hic2.fq,diploid,none\n', '', ['short'], 0, 0)
    assert shortHic.samples[0][0].assembler == 'spades'
    assert shortHic.capabilities.hic && !shortHic.capabilities.scaffold
    def bothReads = check(header + 'both,hifi.bam,one.fq,two.fq,,,diploid,none\n', '', ['both'], 0, 0)
    assert bothReads.samples[0][0].assembler == 'hifiasm'
    check(header + 'good,hifi.bam,,,,,haploid,none\nbad,,one.fq,,,,diploid,none\n', '', ['good'], 1, 0)
    check(header + 'good,hifi.bam,,,,,haploid,none\nbad,,missing.fq,two.fq,,,diploid,none\n', '', ['good'], 1, 0)
    // Taxa stay attached to their own samples; no cohort-level inheritance.
    def species = check('sample_id,hifi_bam,sr_r1,sr_r2,taxid,ploidy\na,hifi.bam,,,1,1\nb,,one.fq,two.fq,2,2\n', '', ['a','b'], 0, 0)
    assert species.samples*.getAt(0)*.taxid == ['1','2']
    log.info 'PASS: actual input parser: mixed routes, singleton, ploidy, Hi-C-only skip, one/two libraries, missing mates/files and species identities'

    // Exercise the real reference process (not its stub), including helper staging.
    def fa = file("${root}/final.fa")
    fa.text = '>chr1\nACGT\n'
    def fai = file("${root}/final.fa.fai")
    fai.text = 'chr1\t4\t6\t4\t5\n'
    def selection = [taxid:'1', species:'Example species', id:'a_primary', sample:'a', haplotype:'primary',
        reference_contigs:['chr1'], fasta:fa.toString(), fai:fai.toString(), selection:'explicit_override', ranking_basis:'fixture']
    SHARED_REFERENCE(Channel.of(tuple('1',selection,fa,fai)),
        file('__REPO__/py_scripts/shared_reference.py'), file('__REPO__/py_scripts/pangenome_input_audit.py'))
    SHARED_REFERENCE.out.manifest.view { taxid, manifest ->
        def result = new groovy.json.JsonSlurper().parse(manifest)
        assert result.bases == 4 && result.coordinate_stage == 'final'
        assert result.sha256.size() == 64
        'PASS: production shared reference process: final FASTA/FAI, checksum and provenance'
    }
}
