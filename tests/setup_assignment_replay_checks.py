"""Tool doubles validate production orchestration, not Cactus assignment biology."""
import gzip
from pathlib import Path
import sys

root=Path(sys.argv[1]);assets=root/'assignment-fixtures';assets.mkdir(exist_ok=True)
tools=root/'assignment-tools';tools.mkdir(exist_ok=True)
packages=tools/'packages';packages.mkdir(exist_ok=True)
for name in ('cactus','cactus/shared','cactus/refmap'):
    directory=packages/name;directory.mkdir(parents=True,exist_ok=True)
    (directory/'__init__.py').write_text('')
metadata=packages/'cactus-3.2.1.dist-info';metadata.mkdir(exist_ok=True)
(metadata/'METADATA').write_text('Name: cactus\nVersion: 3.2.1\n')
(packages/'cactus/shared/version.py').write_text("cactus_commit='d7e5fc1b51327a5a4cdb494c85b0607e3d81a713'\n")
(packages/'cactus/shared/configWrapper.py').write_text('class ConfigWrapper:\n def __init__(self,root):self.xmlRoot=root\n')
(packages/'cactus/shared/common.py').write_text('import os\ndef setupBinaries(options):\n assert options.binariesMode=="local"\n os.environ["CACTUS_BINARIES_MODE"]="local"\n')
(packages/'cactus/refmap/cactus_graphmap.py').write_text('import os\ndef filter_paf(job,source,config,reference=None):\n assert os.environ.get("CACTUS_BINARIES_MODE")=="local", "Cactus library would attempt nested Docker"\n return source\n')
(packages/'cactus/refmap/cactus_minigraph.py').write_text('def minigraph_gfa_from_pansn(job,names,path,source):\n return source\n')
(packages/'cactus/cactus_progressive_config.xml').write_text('''<config><graphmap minMAPQ="5"/><graphmap_split minQueryCoverages="0.75 0.5 0.25" minQueryCoverageThresholds="100000 1000000" minQueryUniqueness="2" ambiguousName="_AMBIGUOUS_"/></config>''')
with gzip.open(assets/'sv.gfa.gz','wt') as f:f.write('S\ts1\tA\n')
(assets/'input.paf').write_text('id=IND_sample.1|composite\t8000000\t0\t4000000\t+\tid=_MINIGRAPH_|s1\t4000000\t0\t4000000\t4000000\t4000000\t60\n')
with gzip.open(assets/'input.gaf.gz','wt') as f:
    f.write((assets/'input.paf').read_text())
(assets/'identities.tsv').write_text('id\tgraph_name\tis_reference\nsample\tIND_sample.1\tFalse\nreference\tREF\tTrue\n')
(assets/'sources.tsv').write_text('assembly\tsequence\tsource_bp\tchromosome_scale\tstatus\nsample\tcomposite\t8000000\tyes\tmapped_but_absent_from_full_and_clip\n')
binary=tools/'rgfa-split'
binary.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
args=sys.argv[1:]
prefix=Path(args[args.index('-b')+1]);log=Path(args[args.index('-L')+1])
regional='-u' in args
chrom='chrA' if regional or args[args.index('-Q')+1]=='1' else '_AMBIGUOUS_'
query='id=IND_sample.1|composite'
log.write_text('tool-double assignment log\\n')
if regional:query+=':1-4000000'
Path(str(prefix)+chrom+'.fa_contigs').write_text(query+'\\n')
Path(str(prefix)+chrom+'.paf').write_text(query+'\\t'+('4000000' if regional else '8000000')+'\\t0\\t4000000\\t+\\ts1\\t4000000\\t0\\t4000000\\t4000000\\t4000000\\t60\\n')
''')
binary.chmod(0o755)
(root/'assignment-check.config').write_text(f'''conda.enabled=false
singularity.enabled=false
docker.enabled=false
trace.enabled=false
report.enabled=false
timeline.enabled=false
dag.enabled=false
params.outdir='{root}/assignment-results'
params.publish_dir_mode='copy'
process {{
 executor='local'
 withName: '.*' {{ cpus=1; memory='1GB'; scratch=false; errorStrategy='terminate'
  beforeScript='export PATH="{tools}:$PATH"; export PYTHONPATH="{packages}"'
 }}
}}
''')
