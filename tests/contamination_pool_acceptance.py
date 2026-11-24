#!/usr/bin/env python3
"""Live Slurm pool test with tiny planted queries and an intentionally queued tail."""
import argparse
from argparse import Namespace
import json
from pathlib import Path
import sys

from contamination_node_fixture import create
from contamination_workers import Results, prepare, digest, atomic_json
from contamination_slurm_pool import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--node', required=True)
    parser.add_argument('--worker-helper', type=Path, required=True)
    args = parser.parse_args()
    if ',' in args.node or '[' in args.node:
        parser.error('Choose exactly one eligible node to exercise a queued second request')
    root = args.root.resolve()
    if not (root / 'plan.json').exists():
        create(root)
        prepare(Namespace(specs=root / 'specifications.json', database=root / 'database/nt',
                          state=root / 'state', threads=1, out=root / 'plan.json'))
    # Keep the successful first worker allocated briefly after publishing results.
    # This leaves the second exclusive request pending while the coordinator retires it.
    wrapper = root / 'linger_worker.py'
    wrapper.write_text('import subprocess,sys,time\n' +
        f'code=subprocess.call([sys.executable,{str(args.worker_helper.resolve())!r},*sys.argv[1:]])\n' +
        'if code==0: time.sleep(20)\nraise SystemExit(code)\n')
    common = dict(plan=root / 'plan.json', helper=wrapper, nodes=2, slots=1,
                  scratch='/scratch', partition='normal', nodelist=args.node, walltime='00:20:00',
                  timeout=1200, poll_seconds=1)
    first = root / 'first/pool_summary.json'
    if not first.exists():
        run(Namespace(out=first.parent, **common))
    first_data = json.loads(first.read_text())
    if len(first_data['jobs']) != 2 or len(first_data['workers']) != 1:
        raise ValueError('Expected one running worker and one retired queued request')
    plan = json.loads((root / 'plan.json').read_text())
    store = Results(plan)
    checksums = {}
    for task in plan['tasks']:
        hits = store.completed(task)
        if not hits:
            raise ValueError('Unverified planted-query result: ' + task['name'])
        index = int(task['name'].split('_')[1].split('.')[0]) - 1
        reference = 'R' + str(index + (12 if task['assembly'] == 'fixture_B' else 0))
        if not any(len(row) == 9 and row[1] == reference and row[2] == '373251' and
                   abs(int(row[4])-int(row[3]))+1 >= 4900 and float(row[7]) >= 99.9
                   for row in (line.split('\t',8) for line in hits.read_text().splitlines())):
            raise ValueError('Expected planted match missing: ' + task['name'])
        checksums[task['key']] = digest(hits)
    run(Namespace(out=root / 'resume', **common))
    resumed = json.loads((root / 'resume/pool_summary.json').read_text())
    if resumed['jobs'] or resumed['remaining']:
        raise ValueError('Complete-cache resume should submit no Slurm jobs')
    atomic_json(root / 'acceptance.json', dict(status='PASS', chunks=len(checksums),
        queued_tail_retired=True, complete_cache_submissions=0, hits_sha256=checksums))


if __name__ == '__main__':
    main()
