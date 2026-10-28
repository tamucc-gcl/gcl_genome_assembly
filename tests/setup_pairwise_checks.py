"""Controlled alignment/QUAST/tidk executables; real evidence and plotting helpers run."""
from pathlib import Path
import sys
import json

root=Path(sys.argv[1]);case=sys.argv[2];repo=Path(__file__).resolve().parents[1]
tools=root/'comparison-tools';tools.mkdir(exist_ok=True)
def tool(name,body):
    p=tools/name;p.write_text('#!'+sys.executable+'\n'+body);p.chmod(0o755)
tool('minimap2','''import sys,os
if '--version' in sys.argv:print('2.30-fixture');sys.exit(0)
case=os.environ['PAIRWISE_CASE']
if case=='alignment_failure':sys.exit(17)
if case=='malformed_paf':print('broken PAF');sys.exit(0)
if case=='no_alignments':sys.exit(0)
print('chr1\\t10000\\t0\\t5000\\t+\\tchr1\\t10000\\t0\\t5000\\t4900\\t5000\\t60')
print('chr1\\t10000\\t6000\\t9000\\t-\\tchr1\\t10000\\t6000\\t9000\\t2940\\t3000\\t60\\tdv:f:0.02')
''')
if case=='native_final':
    folder=root/'native-comparison-fixtures';folder.mkdir(exist_ok=True)
    rows=[]
    for identity,sample,taxid in [('left','specimen','1'),('right','specimen','1'),('singleton','other','2')]:
        fa=folder/(identity+'.fasta');fa.write_text('>chr1\n'+'A'*10000+'\n')
        fai=folder/(identity+'.fasta.fai');fai.write_text('chr1\t10000\t6\t10000\t10001\n')
        rows.append(dict(id=identity,sample=sample,taxid=taxid,completion='finalized',fasta=str(fa),fai=str(fai)))
    manifest=folder/'assembly_eligibility.json';manifest.write_text(json.dumps(dict(assemblies=rows)))
tool('quast.py','''import sys
from pathlib import Path
if '--version' in sys.argv:print('QUAST vfixture');sys.exit(0)
assemblies=sys.argv[1:sys.argv.index('--threads')]
labels=sys.argv[sys.argv.index('--labels')+1].split(',')
assert len(assemblies)==len(labels)>0
assert all(Path(p).is_file() for p in assemblies)
out=Path(sys.argv[sys.argv.index('--output-dir')+1]);out.mkdir()
for suffix in ['tsv','html','txt']:(out/('report.'+suffix)).write_text('fixture')
''')
tool('tidk','''import sys
from pathlib import Path
if '--version' in sys.argv:print('tidk fixture');sys.exit(0)
if sys.argv[1]=='explore':print('id\\trepeat\\tcount\\nchr1\\tCCCTAA\\t1')
elif sys.argv[1]=='search':
 prefix=sys.argv[sys.argv.index('--output')+1];ext=sys.argv[sys.argv.index('--extension')+1]
 Path(prefix+'_telomeric_repeat_windows.'+ext).write_text('id\\twindow\\ttelomeric_repeat_number\\nchr1\\t0\\t5\\n')
elif sys.argv[1]=='plot':Path(sys.argv[sys.argv.index('--output')+1]+'.svg').write_text('<svg/>')
''')
deps=repo/'scratch/comparison-linux-deps'
(root/'pairwise-check.config').write_text(f'''
conda.enabled=false
singularity.enabled=false
docker.enabled=false
trace.enabled=false
report.enabled=false
timeline.enabled=false
dag.enabled=false
params.outdir='{root}/comparison-results/{case}'
params.publish_dir_mode='copy'
params.pairwise_alignment_mode='{'all' if case=='all_species' else ('invalid' if case=='invalid_mode' else 'within_sample')}'
params.run_pairwise_alignments={str(case!='disabled').lower()}
params.pairwise_alignment_min_aln_bp=1
params.pairwise_dotplot_min_seq_bp=0
params.riparian_min_seq_bp=0
params.riparian_min_aln_bp=1
params.pairwise_dotplot_width=5
params.pairwise_dotplot_height=4
params.riparian_width=6
params.riparian_height=4
process {{
 executor='local'
 withName: '.*' {{ cpus=1; memory='1GB'; scratch=false
  beforeScript='export PATH="{tools}:$PATH"; export PYTHONPATH="{deps}"; export MPLCONFIGDIR="{root}/mpl-cache"; export PAIRWISE_CASE="{case}"'
 }}
}}
''')

if case=='native_final':
    with (root/'pairwise-check.config').open('a') as config:config.write(f"\nparams.assembly_manifest='{manifest}'\n")
