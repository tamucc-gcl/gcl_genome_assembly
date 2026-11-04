"""Fixture tools validate orchestration; real R rendering and painting projector run."""
import json
from pathlib import Path
import sys

root=Path(sys.argv[1]);case=sys.argv[2]
fixtures=root/'main-report-fixtures';fixtures.mkdir(exist_ok=True)
tools=root/'main-report-tools';tools.mkdir(exist_ok=True)
minimap=tools/'minimap2'
minimap.write_text('''#!/usr/bin/env python3
import sys
if '--version' in sys.argv:print('2.30-fixture');raise SystemExit
assert '--secondary=no' in sys.argv and '-c' not in sys.argv
assert ('sample_B_' in sys.argv[-1])==('sample_B_' in sys.argv[-2]), 'cross-species painting input'
print('chr1_1\\t1000\\t0\\t1000\\t+\\tchr1_1\\t1000\\t0\\t1000\\t1000\\t1000\\t60\\ttp:A:P')
''');minimap.chmod(0o755)
assemblies=[]
for hap in ['hap1','hap2']:
    aid='sample_'+hap;fa=fixtures/(aid+'.fasta');fai=fixtures/(aid+'.fasta.fai')
    fa.write_text('>chr1_1\n'+'A'*1000+'\n');fai.write_text('chr1_1\t1000\t8\t1000\t1001\n')
    (fixtures/(aid+'.name_map.tsv')).write_text('old_name\tnew_name\tlength\tclass\nold\tchr1_1\t1000\tscaffold\n')
    assemblies.append(dict(id=aid,sample='sample',taxid='1',haplotype=hap,completion='finalized',fasta=str(fa),fai=str(fai),graph_member=False))
references=[assemblies[0]]
if case=='multispecies':
    second=[]
    for source in assemblies:
        row=dict(source);row.update(id='sample_B_'+source['haplotype'],sample='sample_B',taxid='2')
        fa=fixtures/(row['id']+'.fasta');fai=fixtures/(row['id']+'.fasta.fai')
        fa.write_bytes(Path(source['fasta']).read_bytes());fai.write_bytes(Path(source['fai']).read_bytes())
        row.update(fasta=str(fa),fai=str(fai));second.append(row)
        (fixtures/(row['id']+'.name_map.tsv')).write_text('old_name\tnew_name\tlength\tclass\nold\tchr1_1\t1000\tscaffold\n')
    references.append(second[0]);assemblies+=second
(fixtures/'eligibility.json').write_text(json.dumps(dict(assemblies=assemblies,linear_references=references)))
(fixtures/'qc.csv').write_text('sample_id,stage,analysis,metric,hap1,hap2,both\nsample,final,quast,N50,1000,1000,NA\n')
(fixtures/'reads.bam').write_text('input validation fixture\n')
(fixtures/'samples.csv').write_text('sample_id,hifi_bam,hic_r1,hic_r2,sr_r1,sr_r2\nsample,'+str(fixtures/'reads.bam')+',,,,\n')
repo=Path(__file__).resolve().parents[1]
(root/'main-report.config').write_text(f'''conda.enabled=false
singularity.enabled=false
docker.enabled=false
trace.enabled=false
report.enabled=false
timeline.enabled=false
dag.enabled=false
params.run_report=true
params.report_chromosome_painting={'false' if case=='disabled' else 'true'}
params.report_top_scaffolds=20
params.qc_mode='final_only'
params.sample_sheet='{fixtures}/samples.csv'
params.hic_readsets=null
params.taxid=1
params.chimera_detect=false
params.pairwise_alignment_preset='asm5'
params.pairwise_alignment_min_mapq=5
params.pairwise_alignment_min_aln_bp=1
params.outdir='{root}/main-report-results/{case}'
params.publish_dir_mode='copy'
process {{ executor='local'; withName: '.*' {{ cpus=1; memory='2 GB'; scratch=false; errorStrategy='terminate'; beforeScript='export PATH={tools}:/tmp/gcl-summary-check/env/bin:$PATH' }} }}
''')
