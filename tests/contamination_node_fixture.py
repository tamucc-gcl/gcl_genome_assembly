#!/usr/bin/env python3
"""Bounded, species-independent Crest acceptance fixture for the node search workers."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from contamination_workers import Results, atomic_json, digest


def create(root):
    root.mkdir(parents=True, exist_ok=True)
    rng = random.Random(1729)
    sequences = [''.join(rng.choices('ACGT', k=6000)) for _ in range(24)]
    reference = root / 'reference.fasta'
    reference.write_text(''.join(f'>R{i}\n{seq}\n' for i, seq in enumerate(sequences)))
    database = root / 'database'
    database.mkdir()
    subprocess.run(['makeblastdb', '-in', str(reference), '-dbtype', 'nucl', '-parse_seqids',
                    '-taxid', '373251', '-out', str(database / 'nt')], check=True)
    specifications = []
    for i, assembly in enumerate(('fixture_A', 'fixture_B')):
        queries = root / assembly
        queries.mkdir()
        (queries / 'source.json').write_text(json.dumps(dict(assembly=assembly)))
        chunks = ['chunk\tqueries\tquery_bp']
        for j in range(12):
            name = f'protein_{j + 1:05d}.fasta'
            (queries / name).write_text(f'>Q{j + 1:09d}\n{sequences[i * 12 + j][500:5500]}\n')
            chunks.append(f'{name}\t1\t5000')
        (queries / 'chunks.tsv').write_text('\n'.join(chunks) + '\n')
        specifications.append(dict(assembly=assembly, queries=str(queries)))
    atomic_json(root / 'specifications.json', specifications)
    atomic_json(root / 'lock_count.json', 0)


def lock_test(root):
    # Two Slurm ranks contend on the same inode on the actual shared filesystem.
    for _ in range(40):
        with (root / 'filesystem_probe.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            count = json.loads((root / 'lock_count.json').read_text())
            time.sleep(.005)
            atomic_json(root / 'lock_count.json', count + 1)


def verify(args):
    plan = json.loads((args.root / 'plan.json').read_text())
    store = Results(plan)
    if json.loads((args.root / 'lock_count.json').read_text()) != 80:
        raise ValueError('Shared filesystem lock test failed')
    checksums = {}
    for task in plan['tasks']:
        hits = store.completed(task)
        if not hits or not hits.stat().st_size:
            raise ValueError('Fixture search missing its valid planted match: ' + task['name'])
        index = int(task['name'].split('_')[1].split('.')[0]) - 1
        reference = 'R' + str(index + (12 if task['assembly'] == 'fixture_B' else 0))
        rows = [line.split('\t', 8) for line in hits.read_text().splitlines()]
        if not any(len(row) == 9 and row[1] == reference and row[2] == '373251'
                   and abs(int(row[4]) - int(row[3])) + 1 >= 4900 and float(row[7]) >= 99.9
                   for row in rows):
            raise ValueError('Search did not recover the expected planted alignment: ' + task['name'])
        checksums[task['key']] = digest(hits)
    summaries = [json.loads(path.read_text()) for path in args.summaries.glob('worker*.json')]
    if len(summaries) != 2 or sum(row['searched'] for row in summaries) != args.expected:
        raise ValueError('Worker search counts differ from the expected exactly-once execution')
    if args.expected == 24:
        atomic_json(args.root / 'first_hits.json', checksums)
        # Remove exactly one completion marker to exercise recovery, retaining its previous hits.
        (store.directory(plan['tasks'][0]) / 'completion.json').unlink()
    elif checksums != json.loads((args.root / 'first_hits.json').read_text()):
        raise ValueError('Fixture hit content changed after recovery')
    atomic_json(args.root / f'acceptance-{args.expected}.json',
                dict(chunks=len(checksums), searches=sum(row['searched'] for row in summaries),
                     hosts=[row.get('completed', [{}])[0].get('hostname') if row.get('completed') else 'reused'
                            for row in summaries], lock_probe='passed', chunk_recovery='passed'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['create', 'lock-test', 'verify'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--summaries', type=Path)
    parser.add_argument('--expected', type=int)
    args = parser.parse_args()
    args.root = args.root.resolve()
    if args.command == 'create':
        create(args.root)
    elif args.command == 'lock-test':
        lock_test(args.root)
    else:
        verify(args)


if __name__ == '__main__':
    main()
