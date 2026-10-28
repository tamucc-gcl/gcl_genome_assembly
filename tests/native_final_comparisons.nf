nextflow.enable.dsl=2
params.assembly_manifest=null
params.comparison_sample=null
include { PAIRWISE_COMPARISONS } from './workflows/pairwise_comparisons.nf'

process NATIVE_COMPARISON_EVIDENCE {
    cpus 1
    memory '1 GB'
    time '10m'
    conda 'conda-forge::python=3.11'
    errorStrategy 'terminate'
    publishDir "${params.outdir}", mode: params.publish_dir_mode
    input:
    path manifest
    path plan
    path inventory
    path summary
    path plot_status, stageAs:'plot_status/*'
    path retention_files, stageAs:'retention/*'
    val selected_ids
    output:
    path 'comparison_acceptance.md'
    script:
    """
    python3 - <<'PY'
import csv,hashlib,json
from pathlib import Path
def table(path):
    with open(path) as handle:return list(csv.DictReader(handle,delimiter='\\t'))
selected=${groovy.json.JsonOutput.toJson(selected_ids)}
members={r['id']:r for r in json.loads(Path('${manifest}').read_text())['assemblies'] if r['id'] in selected and r['completion']=='finalized'}
inputs=table('${inventory}');assert len(inputs)==len(members)
assert {r['id'] for r in inputs}==set(members)
for row in inputs:
    source=members[row['id']]
    assert row['coordinate_stage']=='final' and row['fasta']==source['fasta']
    assert hashlib.sha256(Path(source['fai']).read_bytes()).hexdigest()==row['fai_sha256']
planned=table('${plan}');measurements=table('${summary}')
keys=lambda rows:{(r['ref_id'],r['query_id'] if 'query_id' in r else r['qry_id']) for r in rows}
assert keys(planned)==keys(measurements) and len(planned)==len(measurements)
for row in planned:
    left,right=members[row['ref_id']],members[row['query_id']]
    assert left['taxid']==right['taxid'] and row['coordinate_stage']=='final'
    if row['mode']=='within_sample':assert left['sample']==right['sample']
retention={}
for path in Path('retention').glob('*.tsv'):
    for row in table(path):
        assert 0<=float(row['ref_coverage_fraction'])<=1 and 0<=float(row['query_coverage_fraction'])<=1
        retention[(row['ref_id'],row['qry_id'],row['filter'])]=row
statuses=[table(p) for p in Path('plot_status').glob('*.tsv')]
assert len(statuses)==len(planned)
assert all({r['plot'] for r in rows}=={'dotplot','ribbon'} and len(rows)==2 for rows in statuses)
text=['# Native final-assembly comparison check','',
      'PASS execution: expected comparisons, evidence summaries and both R plot types completed.',
      f'{len(inputs)} finalized assemblies; {len(planned)} comparisons. Final index fingerprints match the source manifest.',
      '','| Reference | Query | Chains | Reference coverage | Query coverage | Approx. match fraction | Figures |','|---|---|---:|---:|---:|---:|---|']
for row in measurements:
    prefix=row['ref_id']+'_vs_'+row['qry_id']
    fraction=row['approx_match_fraction']
    assert row['measurement_basis']=='mapping_chain_approximate'
    cov=retention[(row['ref_id'],row['qry_id'],'min_'+row['min_aln_bp'])]
    assert cov['n_align']==row['n_align']
    coverage=f"{float(cov['ref_coverage_fraction']):.1%} | {float(cov['query_coverage_fraction']):.1%}"
    text.append('| '+row['ref_id']+' | '+row['qry_id']+' | '+row['n_align']+' | '+coverage+' | '+fraction+' | '+
                '[Dotplot](pairwise_alignments/'+prefix+'_dotplot.png) / '+
                '[Ribbons](pairwise_alignments/'+prefix+'_riparian.png) |')
text.extend(['','Mapping-chain match counts and block lengths are approximate; their ratio is not measured base-level or whole-genome identity.',
             'Coverage is the union of mapped coordinate intervals divided by full final assembly length. It is chain coverage, not base-by-base alignment.',
             'Per-pair retention tables show unfiltered, MAPQ-only and length-filtered coverage; raw mappings remain on Crest.',
             'Execution PASS does not assert that mapping coverage is biologically sufficient.',
             'Axes include every final sequence meeting the plotting size threshold, including unaligned sequences.',
             'Empty plots have explicit status tables. Biological assembly quality is not certified by this execution check.'])
Path('comparison_acceptance.md').write_text('\\n'.join(text)+'\\n')
PY
    """
}

workflow {
    if (!params.assembly_manifest) error 'Supply --assembly_manifest from the accepted finalized assembly run'
    def sourceManifest=file(params.assembly_manifest,checkIfExists:true)
    def data=new groovy.json.JsonSlurper().parse(sourceManifest)
    def members=data.assemblies.findAll { it.completion=='finalized' }
    if (params.comparison_sample) members=members.findAll { it.sample==params.comparison_sample }
    if (!members) error 'No finalized assemblies in the supplied manifest'
    records=members.collect { row ->
        tuple([id:row.id,sample:row.sample,taxid:row.taxid],
              file(row.fasta,checkIfExists:true),file(row.fai,checkIfExists:true))
    }
    PAIRWISE_COMPARISONS(Channel.fromList(records))
    NATIVE_COMPARISON_EVIDENCE(sourceManifest,PAIRWISE_COMPARISONS.out.plan,
        PAIRWISE_COMPARISONS.out.inventory,PAIRWISE_COMPARISONS.out.summary,
        PAIRWISE_COMPARISONS.out.status.map { a,b,status -> status }.toList(),
        PAIRWISE_COMPARISONS.out.retention.map { a,b,retention -> retention }.toList(),members.collect { it.id })
}
