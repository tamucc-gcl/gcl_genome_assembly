"""External CLI doubles; real workflows, validators and GenBank converter execute."""
from pathlib import Path
import os
import sys

root=Path(sys.argv[1]);case=sys.argv[2];repo=Path(__file__).resolve().parents[1]
tools=root/'organelle-tools';tools.mkdir(exist_ok=True)
real=sys.executable
def tool(name,body):
    p=tools/name;p.write_text('#!'+real+'\n'+body);p.chmod(0o755)

body='''import sys,os
from pathlib import Path
case=os.environ['ORG_CASE']
if '--version' in sys.argv:print('GetOrganelle v1.7.7.1');sys.exit(0)
out=Path(sys.argv[sys.argv.index('-o')+1]);out.mkdir(exist_ok=True)
org=sys.argv[sys.argv.index('-F')+1]
retry='from_assembly' in Path(sys.argv[0]).name
if case=='execution_failure':print('native execution failed');sys.exit(17)
if case=='logged_error':print('2026-10-25 - ERROR: missing dependency');print('Thank you!');sys.exit(0)
if case=='incomplete_log':print('Interrupted before completion');sys.exit(0)
if case in ('retry','fungus_retry') and not retry:
    graph=out/'extended_spades/K55/assembly_graph.fastg';graph.parent.mkdir(parents=True,exist_ok=True);graph.write_text('graph\\n')
elif case!='no_recovery':
    status='incomplete' if case=='incomplete' else ('scaffolds' if case in ('linear','annotation_failure','missing_bed','multiple_records') else 'complete')
    fasta='>contig\\nACGTACGTACGTACGT\\n'
    if case=='multiple_records':fasta+='>second\\nACGT\\n'
    (out/(org+'.K55.'+status+'.graph1.1.path_sequence.fasta')).write_text(fasta)
(out/'get_org.log.txt').write_text('Completed fixture\\n')
print('Thank you!')
'''
tool('get_organelle_from_reads.py',body);tool('get_organelle_from_assembly.py',body)
tool('get_organelle_config.py',"import sys\nprint('download failed')\nsys.exit(17)\n")
tool('curl',"import sys\nprint('fixture database download failure')\nsys.exit(17)\n")
tool('runmitos','''import sys,os
from pathlib import Path
if '--version' in sys.argv:print('MITOS2 fixture');sys.exit(0)
if os.environ['ORG_CASE']=='annotation_failure':sys.exit(19)
out=Path(sys.argv[sys.argv.index('-o')+1]);out.mkdir(exist_ok=True)
if os.environ['ORG_CASE']!='missing_bed':(out/'result.bed').write_text('contig\\t0\\t9\\tcox1\\t0\\t+\\n')
''')
tool('taxonkit','''import sys,os
mode=os.environ['ORG_CASE'];taxid=sys.stdin.read().strip()
if mode=='taxonomy_failure':sys.exit(17)
if mode=='unknown_taxid':print(taxid);sys.exit(0)
if sys.argv[1]=='reformat':print(taxid+'\\tEukaryota\\tMetazoa\\tChordata\\tActinopteri\\tOrder\\tFamily\\tGenus\\tSpecies')
elif mode!='empty_lineage':print(taxid+'\\tcellular organisms;Eukaryota;Metazoa;Chordata;Actinopteri')
''')
# Plot rendering is an external boundary double; conversion runs with real Biopython.
python_body='''import sys,os
from pathlib import Path
if len(sys.argv)>1 and Path(sys.argv[1]).name=='plot_mito_circular.py':
    Path(sys.argv[sys.argv.index('--output')+1]).write_bytes(b'PNG fixture')
else:os.execv(REAL,[REAL]+sys.argv[1:])
'''.replace('REAL',repr(real))
tool('python3',python_body);tool('python',python_body)
db=root/('go-db-'+case)
if case!='database_failure':
    for folder in ('SeedDatabase','LabelDatabase'):
        p=db/folder;p.mkdir(parents=True,exist_ok=True);(p/'VERSION').write_text('fixture\n')
        for name in ('embplant_pt','embplant_mt','embplant_nr','fungus_mt','fungus_nr','animal_mt','other_pt'):
            (p/(name+'.fasta')).write_text('>seed\nACGT\n')
mitos_parent=root/('mitos-broken' if case=='mitos_database_failure' else 'mitos-fixture')
mitos=mitos_parent/'refseq89m'
if case!='mitos_database_failure':
    mitos.mkdir(parents=True,exist_ok=True);(mitos/'data').write_text('fixture')
    fungus=mitos_parent/'refseq89f';fungus.mkdir(exist_ok=True);(fungus/'data').write_text('fixture')
deps=repo/'scratch/organelle-test-deps'
(root/'organelle-check.config').write_text('''
conda.enabled=false
singularity.enabled=false
docker.enabled=false
trace.enabled=false
report.enabled=false
timeline.enabled=false
dag.enabled=false
params.outdir='%s/organelle-results/%s'
params.publish_dir_mode='copy'
params.getorganelle_downloads='%s'
params.mitos_downloads='%s'
process {
    executor='local'
    withName: '.*' { cpus=1; memory='1GB'; scratch=false; errorStrategy='terminate'
        beforeScript='export PATH="%s:$PATH"; export PYTHONPATH="%s"; export ORG_CASE="%s"'
    }
}
'''%(root,case,db,mitos_parent,tools,deps,case))
