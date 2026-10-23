"""Tiny tool doubles exercise the production diagnostic workflow, not real Cactus biology."""
import gzip
from pathlib import Path
import runpy
import sys

root=Path(sys.argv[1])
runpy.run_path(str(Path(__file__).with_name('setup_assignment_replay_checks.py')),run_name='__main__')
assets=root/'deletion-fixtures';assets.mkdir(exist_ok=True)
tools=root/'assignment-tools';packages=tools/'packages'
(packages/'cactus/cactus_progressive_config.xml').write_text('''<config><graphmap minMAPQ="5" collapse="none" delFilter="10000000" delFilterThreshold="0.01" delFilterQuerySizeThreshold="2"/><graphmap_split minQueryCoverages="0.75 0.5 0.25" minQueryCoverageThresholds="100000 1000000" minQueryUniqueness="2" ambiguousName="_AMBIGUOUS_"/></config>''')
with gzip.open(assets/'sv.gfa.gz','wt') as f:
    f.write('S\ts1\tACGT\tSN:Z:chrA\tSO:i:0\tSR:i:0\nS\ts2\tACGT\tSN:Z:chrB\tSO:i:0\tSR:i:0\n')
query='id=IND_sample.1|chrA+chrB'
records=[f'{query}\t8000000\t{start}\t{start+4000000}\t+\tid=_MINIGRAPH_|{node}\t4000000\t0\t4000000\t4000000\t4000000\t60\n'
         for start,node in ((0,'s1'),(4000000,'s2'))]
with gzip.open(assets/'input.paf.unfiltered.gz','wt') as f:f.writelines(records)
(assets/'input.paf').write_text(records[1])
(assets/'input.paf.filter.log').write_text('fixture: native drops the chrA lane from the composite\n')
(assets/'identities.tsv').write_text('id\tgraph_name\tis_reference\nsample\tIND_sample.1\tFalse\nreference\tREF\tTrue\n')
(assets/'sources.tsv').write_text('assembly\tsequence\tsource_bp\tchromosome_scale\tstatus\tclip_retained_path_bp\nsample\tchrA+chrB\t8000000\tyes\trepresented_in_clip\t4000000\n')
binary=tools/'vg'
binary.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
args=sys.argv[1:]
Path(args[args.index('-T')+1]).write_text('fixture translation\\n')
print('fixture graph')
''');binary.chmod(0o755)
binary=tools/'filter-paf-deletions'
binary.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
from collections import defaultdict
args=sys.argv[1:];rows=Path(args[2]).read_text().splitlines(True)
targets=defaultdict(set)
for line in rows:
 f=line.split('\\t');targets[f[0]].add(f[5])
for line in rows:
 f=line.split('\\t')
 if len(targets[f[0]])>1 and f[5].endswith('s1'):continue
 sys.stdout.write(line)
print('fixture native deletion filter',file=sys.stderr)
''');binary.chmod(0o755)
binary=tools/'rgfa-split'
binary.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
args=sys.argv[1:];prefix=args[args.index('-b')+1]
rows=[line.rstrip().split('\\t') for line in Path(args[args.index('-p')+1]).read_text().splitlines()]
regional='-u' in args
Path(args[args.index('-L')+1]).write_text('fixture assignment decisions\\n')
for chrom,node in [('chrA','s1'),('chrB','s2')]:
 chosen=[f for f in rows if f[5].endswith(node)]
 if not regional and chrom!='chrB':continue
 if not chosen:continue
 query=chosen[0][0];lo=min(int(f[2]) for f in chosen);hi=max(int(f[3]) for f in chosen)
 name=query+(':'+str(lo+1)+'-'+str(hi) if regional else '')
 Path(prefix+chrom+'.fa_contigs').write_text(name+'\\n')
 with Path(prefix+chrom+'.paf').open('w') as out:
  for f in chosen:
   f[0]=name
   if regional:f[1]=str(hi-lo);f[2]=str(int(f[2])-lo);f[3]=str(int(f[3])-lo)
   out.write('\\t'.join(f)+'\\n')
''');binary.chmod(0o755)
(root/'deletion-check.config').write_text(f'''conda.enabled=false
singularity.enabled=false
docker.enabled=false
trace.enabled=false
report.enabled=false
timeline.enabled=false
dag.enabled=false
params.outdir='{root}/deletion-results'
params.publish_dir_mode='copy'
process {{
 executor='local'
 withName: '.*' {{ cpus=1; memory='1GB'; scratch=false; errorStrategy='terminate'
  beforeScript='export PATH="{tools}:$PATH"; export PYTHONPATH="{packages}"'
 }}
}}
''')
