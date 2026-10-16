"""Tiny external tool doubles; orchestration/auditors remain production code."""
import sys
from pathlib import Path

root = Path(sys.argv[1])
case = sys.argv[2]
tools = root / 'core-tools'
tools.mkdir(exist_ok=True)
double = tools / 'cactus-pangenome'
double.write_text('''#!/usr/bin/env python3
import gzip,sys
from pathlib import Path
if '--help' in sys.argv:
    print('--gref --haplo --refContigs --vcfbub --restart --clean'); sys.exit(0)
out = Path(sys.argv[sys.argv.index('--outDir')+1]); out.mkdir(exist_ok=True)
taxid = sys.argv[sys.argv.index('--outName')+1]
for suffix in ['gbz','og','hapl','snarls','gref.gbz','stats.tgz']:
    (out/(taxid+'.'+suffix)).write_text('tool double artifact\\n')
for suffix in ['gfa.gz','gref.gfa.gz','full.gfa.gz']:
    with gzip.open(out/(taxid+'.'+suffix),'wt') as f: f.write('H\\tVN:Z:1.1\\nS\\t1\\tACGT\\n')
with gzip.open(out/(taxid+'.gref.vcf.gz'),'wt') as f:
    f.write('##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n')
(out/(taxid+'.gref.vcf.gz.tbi')).write_text('index tool double\\n')
''')
if case == 'missing_graph_role':
    double.write_text(double.read_text() + "\n(out/(taxid+'.gref.gbz')).unlink()\n")
double.chmod(0o755)
if case == 'graph_failure':
    double.write_text(double.read_text().replace("out = Path(sys.argv", "sys.exit(17)\nout = Path(sys.argv"))
cactus = tools / 'cactus'
cactus.write_text('#!/bin/sh\necho "cactus tool double"\n')
cactus.chmod(0o755)
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
        beforeScript='export PATH="%s:$PATH"'
    }
}
''' % (root, case, tools))
