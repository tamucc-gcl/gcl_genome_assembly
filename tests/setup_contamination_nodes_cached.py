"""Prime a tiny synthetic chunk cache for the real no-worker Nextflow branch."""
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

repo = Path(sys.argv[1])
root = Path('/tmp/contamination-nodes-cached')
root.mkdir(exist_ok=True)
sys.path.insert(0, str(repo / 'py_scripts'))
from contamination_workers import prepare, Results, digest, atomic_json

tools = root / 'bin'
tools.mkdir(exist_ok=True)
blast = tools / 'blastn'
blast.write_text('#!/bin/sh\nif [ "$1" = -version ]; then echo "blastn: 2.16.0+"; else exit 99; fi\n')
blast.chmod(0o755)
os.environ['PATH'] = str(tools) + ':' + os.environ['PATH']
database = root / 'db'
database.mkdir(exist_ok=True)
(database / 'nt.nal').write_text('DBLIST nt.000\n')
(database / 'nt.000.ndb').write_text('synthetic database fixture')
specs = []
for assembly in ('A', 'B'):
    queries = root / assembly
    queries.mkdir(exist_ok=True)
    (queries / 'source.json').write_text(json.dumps(dict(assembly=assembly)))
    (queries / 'chunks.tsv').write_text('chunk\tqueries\tquery_bp\nchunk.fasta\t1\t4\n')
    (queries / 'chunk.fasta').write_text('>' + assembly + '\nACGT\n')
    specs.append(dict(assembly=assembly, queries=str(queries)))
(root / 'specs.json').write_text(json.dumps(specs))
prepare(SimpleNamespace(specs=root / 'specs.json', database=database / 'nt',
                        state=Path('/tmp/chimera-nxf-runtime/work/contamination_megablast'),
                        threads=1, out=root / 'plan.json'))
plan = json.loads((root / 'plan.json').read_text())
store = Results(plan)
for task in plan['tasks']:
    directory = store.directory(task)
    directory.mkdir(parents=True, exist_ok=True)
    hits = directory / 'hits.tsv'
    hits.write_text(task['assembly'] + '\tfixture\n')
    atomic_json(directory / 'completion.json', dict(key=task['key'], query_sha256=task['query_sha256'],
                                                  hits_bytes=hits.stat().st_size, hits_sha256=digest(hits)))
