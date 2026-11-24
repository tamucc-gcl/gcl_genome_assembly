#!/usr/bin/env python3
"""Verify completed chunk evidence and preserve measurements before a clean resume."""
import argparse
import json
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'py_scripts'))
from contamination_workers import Results

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--project', type=Path, required=True)
parser.add_argument('--plan', type=Path, required=True)
parser.add_argument('--controller', required=True)
args = parser.parse_args()
project = args.project.resolve(strict=True)
run = project / 'comparisons/contamination-staged-cohort-01'
plan = json.loads(args.plan.read_text())
expected = json.loads((run / 'expected.json').read_text())
if set(plan['assemblies']) != set(expected) or Path(plan['state']).resolve() != (project / 'work/contamination_megablast').resolve():
    raise ValueError('Plan does not belong to this cohort/cache')
store = Results(plan)
missing = [task for task in plan['tasks'] if not store.completed(task)]
print(f"Checksum-verified completed chunks: {len(plan['tasks']) - len(missing)} / {len(plan['tasks'])}", flush=True)
if missing:
    raise SystemExit('Work remains; no cancellation or restart is authorized by this recovery script.')
out = run / 'recovery-before-resume'
out.mkdir(exist_ok=True)
keys = {task['key'] for task in plan['tasks']}
saved = []
for path in (project / 'work').glob('*/*/worker_summary.json'):
    data = json.loads(path.read_text())
    if any(record.get('key') in keys for record in data.get('completed', [])):
        target = out / (path.parent.name + '.worker_summary.json')
        if not target.exists():
            shutil.copyfile(path, target)
        saved.append(target.name)
for path in (args.plan, run / f'nextflow-{args.controller}.log', run / f'slurm-{args.controller}.out',
             run / 'results/pipeline/pipeline_trace.txt'):
    if path.is_file() and not (out / path.name).exists():
        shutil.copyfile(path, out / path.name)
(out / 'verified.json').write_text(json.dumps(dict(chunks=len(plan['tasks']), missing=[], worker_summaries=saved), indent=2) + '\n')
print(f'Preserved {len(saved)} completed-worker summaries and original logs in {out}', flush=True)
