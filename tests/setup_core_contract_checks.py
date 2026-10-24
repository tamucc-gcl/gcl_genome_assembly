"""Tiny external tool doubles; orchestration/auditors remain production code."""
import sys
import shutil
from pathlib import Path

root = Path(sys.argv[1])
case = sys.argv[2]
tools = root / 'core-tools'
tools.mkdir(exist_ok=True)
packages=tools/'packages';(packages/'cactus/shared').mkdir(parents=True,exist_ok=True)
(packages/'cactus/__init__.py').write_text('')
(packages/'cactus/shared/__init__.py').write_text('')
(packages/'cactus/shared/version.py').write_text("cactus_commit='d7e5fc1b51327a5a4cdb494c85b0607e3d81a713'\n")
shutil.copyfile(Path(__file__).parent/'fixtures/cactus_fragment_rename.py',packages/'cactus/shared/common.py')
# External extraction/sanitizer doubles support the native header probe. The
# unchanged upstream sed helper above executes for real; native tools run on Crest.
for name,body in {
    'samtools': '''import sys
from pathlib import Path
fa=Path(sys.argv[2]).read_text().splitlines();name=fa[0][1:];sequence=''.join(fa[1:])
for region in Path(sys.argv[sys.argv.index('--region-file')+1]).read_text().splitlines():
    base,span=region.rsplit(':',1);lo,hi=map(int,span.split('-'))
    assert base==name
    print('>'+region);print(sequence[lo-1:hi])
''',
    'cactus_sanitizeFastaHeaders': '''import re,sys
from pathlib import Path
for line in Path(sys.argv[1]).read_text().splitlines():
    if line.startswith('>'):
        match=re.fullmatch(r'(.*):(\\d+)-(\\d+)',line)
        if match:line=match[1]+'_sub_'+str(int(match[2])-1)+'_'+match[3]
    print(line)
'''
}.items():
    tool=tools/name;tool.write_text('#!/usr/bin/env python3\n'+body);tool.chmod(0o755)
if case=='fragment_header_failure':
    sanitizer=tools/'cactus_sanitizeFastaHeaders'
    sanitizer.write_text(sanitizer.read_text().replace("print(line)","print(line.replace('_sub_', '_wrong_'))"))
(packages/'cactus/cactus_progressive_config.xml').write_text('<config><graphmap collapse="none" delFilter="10000000"/><graphmap_split ambiguousName="_AMBIGUOUS_"/></config>')
(packages/'cactus-3.2.1.dist-info').mkdir(exist_ok=True)
(packages/'cactus-3.2.1.dist-info/METADATA').write_text('Name: cactus\nVersion: 3.2.1\n')
filter_double=tools/'filter-paf-deletions'
filter_double.write_text('#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\nsys.stdout.write(Path(sys.argv[3]).read_text())\n')
filter_double.chmod(0o755)
double = tools / 'cactus-pangenome'
double.write_text('''#!/usr/bin/env python3
import gzip,sys,tarfile,io,subprocess
from pathlib import Path
if '--help' in sys.argv:
    print('--gref --haplo --refContigs --vcfbub --restart --clean'); sys.exit(0)
out = Path(sys.argv[sys.argv.index('--outDir')+1]); out.mkdir(exist_ok=True)
taxid = sys.argv[sys.argv.index('--outName')+1]
# Exercise the production wrappers from a child process, as Cactus worker jobs do.
scratch=Path('fixture-worker');scratch.mkdir()
(scratch/'mg.gfa').write_text('S\\ts1\\tACGT\\tSN:Z:id=REF_'+taxid+'|chr1\\tSO:i:0\\tSR:i:0\\n')
(scratch/'mg.gfa.vg').write_text('fixture graph\\n');(scratch/'trans').write_text('fixture trans\\n')
(scratch/'raw.paf').write_text('id=IND_'+taxid+'|chr1\\t4\\t0\\t4\\t+\\ts1\\t4\\t0\\t4\\t4\\t4\\t60\\n')
with (scratch/'guarded.paf').open('w') as f:
    subprocess.run(['filter-paf-deletions',str(scratch/'mg.gfa.vg'),str(scratch/'trans'),str(scratch/'raw.paf'),'-d','10000000'],stdout=f,check=True)
subprocess.run(['rgfa-split','-g',str(scratch/'mg.gfa'),'-p',str(scratch/'guarded.paf'),'-b',str(scratch/'split_'),
                '-a','_AMBIGUOUS_','-n','.25','-Q','2','-G','-L',str(scratch/'split.log')],check=True)
for suffix in ['gbz','og','hapl','snarls','gref.gbz','stats.tgz']:
    (out/(taxid+'.'+suffix)).write_text('tool double artifact\\n')
for suffix in ['gfa.gz','gref.gfa.gz','full.gfa.gz','sv.gfa.gz']:
    with gzip.open(out/(taxid+'.'+suffix),'wt') as f:
        f.write('H\\tVN:Z:1.1\\nS\\t1\\tACGT\\nW\\tREF_'+taxid+'\\t0\\tchr1\\t0\\t4\\t>1\\nW\\tIND_'+taxid+'\\t0\\tchr1\\t0\\t4\\t>1\\n')
with gzip.open(out/(taxid+'.gref.vcf.gz'),'wt') as f:
    f.write('##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n')
(out/(taxid+'.gref.vcf.gz.tbi')).write_text('index tool double\\n')
mapping='id=IND_'+taxid+'|chr1\\t4\\t0\\t4\\t+\\tid=_MINIGRAPH_|s1\\t4\\t0\\t4\\t4\\t4\\t60\\n'
with gzip.open(out/(taxid+'.gaf.gz'),'wt') as f: f.write(mapping)
(out/(taxid+'.paf')).write_text(mapping)
with gzip.open(out/(taxid+'.paf.unfiltered.gz'),'wt') as f: f.write(mapping)
(out/(taxid+'.paf.filter.log')).write_text('fixture deletion filter\\n')
with tarfile.open(out/(taxid+'.stats.tgz'),'w:gz') as archive:
    text='chr1\\tREF_'+taxid+'#0\\t4\\nchr1\\tIND_'+taxid+'#0\\t4\\n'
    paths='chr1\\tREF_'+taxid+'#0#chr1\\t4\\nchr1\\tIND_'+taxid+'#0#chr1\\t4\\n'
    for name,body in [('sample-stats.tsv',text),('path-stats.tsv',paths),('clip-stats.bed','')]:
        data=body.encode(); member=tarfile.TarInfo(name); member.size=len(data)
        archive.addfile(member,io.BytesIO(data))
''')
if case == 'missing_graph_role':
    double.write_text(double.read_text() + "\n(out/(taxid+'.gref.gbz')).unlink()\n")
double.chmod(0o755)
if case == 'graph_failure':
    double.write_text(double.read_text().replace("out = Path(sys.argv", "sys.exit(17)\nout = Path(sys.argv"))
cactus = tools / 'cactus'
cactus.write_text('#!/bin/sh\necho "cactus tool double"\n')
cactus.chmod(0o755)
for name, body in [('vg','nodes\\t1\\nedges\\t0\\nlength\\t4\\n'),
                   ('odgi','#length\\tnodes\\tedges\\tpaths\\tsteps\\n4\\t1\\t0\\t2\\t2\\n')]:
    tool = tools/name
    tool.write_text("#!/bin/sh\nprintf '"+body+"'\n")
    tool.chmod(0o755)
mito = tools / 'mitohifi.py'
mito.write_text('#!/bin/sh\n[ "$1" != "-v" ] || { echo "mitohifi tool double"; exit 0; }\n' +
                ('exit 17\n' if case == 'organelle_failure' else 'exit 0\n'))
mito.chmod(0o755)
(root / 'core-check.config').write_text('''conda.enabled=false
singularity.enabled=false
docker.enabled=false
trace.enabled=false
report.enabled=false
timeline.enabled=false
dag.enabled=false
params.outdir="%s/core-results"
params.publish_dir_mode='copy'
params.contract_case='%s'
process {
    executor='local'
    withName: '.*' { cpus=1; memory='1GB'; scratch=false; errorStrategy='terminate'
        beforeScript='export PATH="%s:$PATH"; export PYTHONPATH="%s"'
    }
}
''' % (root, case, tools,packages))
with (root/'core-check.config').open('a') as handle:
    handle.write('params.pangenome_assignment_replay='+('true' if case=='assignment_inputs' else 'false')+'\n')
    handle.write('params.pangenome_deletion_filter_replay='+('true' if case=='deletion_inputs' else 'false')+'\n')
