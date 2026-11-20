#!/usr/bin/env python3
"""Benchmark unchanged DIAMOND evidence searches on shared storage and node-local scratch."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import statistics
import subprocess
import tempfile
import time

FLAGS = ['--sensitive', '--evalue', '1e-5', '--frameshift', '15', '--top', '20',
         '--range-culling', '--outfmt', '6', 'qseqid', 'sseqid', 'staxids',
         'qstart', 'qend', 'bitscore', 'evalue', 'pident', 'stitle']


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n')


def select(project, trace, count, output):
    by_assembly = {}
    with trace.open(newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            match = re.search(r'(?:^|:)CONTAMINATION_PREPARE \(([^)]+)\)$', row['name'])
            if not match or row['status'] not in ('CACHED', 'COMPLETED'):
                continue
            prefix, partial = row['hash'].split('/')
            tasks = list((project / 'work' / prefix).glob(partial + '*'))
            if len(tasks) != 1 or not (tasks[0] / 'queries/chunks.tsv').is_file():
                raise ValueError('Cannot resolve retained prepared queries for ' + match[1])
            directory = tasks[0] / 'queries'
            with (directory / 'chunks.tsv').open(newline='') as source:
                chunks = list(csv.DictReader(source, delimiter='\t'))
            # Spread probes across each assembly, including beginning, middle and tail.
            order = sorted(range(len(chunks)), key=lambda i: (min(i, abs(i-len(chunks)//2), len(chunks)-1-i), i))
            by_assembly[match[1]] = [(directory / chunks[i]['chunk'], int(chunks[i]['query_bp'])) for i in order]
    if not by_assembly:
        raise ValueError('No retained CONTAMINATION_PREPARE tasks found in the supplied trace')
    for index, assembly in enumerate(sorted(by_assembly)):
        values = by_assembly[assembly]
        offset = index % min(3, len(values)) if values else 0
        by_assembly[assembly] = values[offset:] + values[:offset]
    candidates = [(assembly, values[i]) for i in range(max(map(len, by_assembly.values())))
                  for assembly, values in sorted(by_assembly.items()) if i < len(values)]
    chosen = candidates[:count]
    if len(chosen) < count:
        raise ValueError(f'Only {len(chosen)} retained chunks available; requested {count}')
    output.mkdir(parents=True)
    rows = []
    for index, (assembly, (source, bp)) in enumerate(chosen):
        name = f'{index:03d}_{assembly}_{source.name}'
        target = output / name
        shutil.copyfile(source, target)
        rows.append(dict(assembly=assembly, file=name, original=str(source.resolve()),
                         query_bp=bp, sha256=sha256(target)))
    write_json(output / 'selection.json', dict(cohort_chunks=len(candidates), samples=rows,
                                              cohort_query_bp=sum(bp for assembly, (path, bp) in candidates),
                                              trace=str(trace), project=str(project)))


def normalized_hits(path):
    return hashlib.sha256(''.join(sorted(path.read_text().splitlines(keepends=True))).encode()).hexdigest()


def benchmark(args):
    selection = json.loads((args.queries / 'selection.json').read_text())
    version = subprocess.check_output(['diamond', 'version'], text=True).strip()
    if version != 'diamond version 2.1.14':
        raise ValueError('Benchmark requires the pipeline-pinned DIAMOND 2.1.14')
    cpus = sorted(os.sched_getaffinity(0))
    if len(cpus) < args.jobs * args.threads:
        raise ValueError('Not enough allocated CPUs for jobs multiplied by threads')
    timing_binary = shutil.which('time')
    if not timing_binary:
        raise ValueError('GNU time is required to record per-search peak memory')
    database = args.bundle / 'reference_proteomes.dmnd'
    identity = json.loads((args.bundle / 'protein_database.json').read_text())
    before = (database.stat().st_size, database.stat().st_mtime_ns)
    needed = before[0] + args.jobs * 100 * 1024**3
    usage = shutil.disk_usage(args.scratch)
    if usage.free < needed:
        raise ValueError(f'{os.uname().nodename}: need {needed/1024**3:.1f} GiB scratch '
                         f'(index plus 100 GiB per search); available {usage.free/1024**3:.1f} GiB')
    args.out.mkdir(parents=True)
    local = Path(tempfile.mkdtemp(prefix='contamination-protein-benchmark-', dir=args.scratch))
    try:
        local_db = local / 'reference_proteomes.dmnd'
        started = time.monotonic()
        shutil.copyfile(database, local_db)
        if sha256(local_db) != identity['database_sha256']:
            raise ValueError('Staged protein database does not match its source manifest')
        staging = time.monotonic() - started
        metrics = []
        for mode, db in [('shared', database), ('scratch', local_db)]:
            print(f'Starting {mode}: {args.jobs} concurrent searches, {args.threads} threads each', flush=True)
            directory = args.out / mode
            directory.mkdir()
            pending = list(selection['samples'])

            def slot(index):
                results = []
                assigned = cpus[index*args.threads:(index+1)*args.threads]
                for sample in pending[index::args.jobs]:
                    query = args.queries / sample['file']
                    if sha256(query) != sample['sha256']:
                        raise ValueError('Benchmark query changed')
                    temporary = local / f'{mode}-{index}'
                    temporary.mkdir(exist_ok=True)
                    hits = directory / (sample['file'] + '.hits.tsv')
                    timings = directory / (sample['file'] + '.time.txt')
                    command = [timing_binary, '-v', '-o', str(timings), 'taskset', '-c',
                               ','.join(map(str, assigned)), 'diamond', 'blastx', '--query', str(query),
                               '--db', str(db), '--out', str(hits), '--threads', str(args.threads),
                               '--tmpdir', str(temporary), *FLAGS]
                    write_json(directory / (sample['file'] + '.command.json'), command)
                    begin = time.monotonic()
                    with (directory / (sample['file'] + '.stderr.txt')).open('w') as stderr:
                        subprocess.run(command, stderr=stderr, check=True)
                    text = timings.read_text()
                    rss = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', text)
                    results.append(dict(mode=mode, file=sample['file'], seconds=time.monotonic()-begin,
                                        peak_rss_kib=int(rss[1]) if rss else None,
                                        normalized_hits_sha256=normalized_hits(hits)))
                    print(f'Completed {mode}/{sample["file"]}', flush=True)
                return results

            phase_start = time.monotonic()
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                observations = [row for rows in pool.map(slot, range(args.jobs)) for row in rows]
            elapsed = time.monotonic() - phase_start
            wall = [row['seconds'] for row in observations]
            ideal = sum(wall)/len(wall)*selection['cohort_chunks']/args.jobs/3600
            per_base = sum(wall) / sum(row['query_bp'] for row in selection['samples'])
            metrics.append(dict(mode=mode, seconds=elapsed, observations=observations,
                                median_chunk_seconds=statistics.median(wall),
                                estimated_cohort_hours_scaled_by_query_bp=per_base*selection['cohort_query_bp']/args.jobs/3600,
                                estimated_cohort_hours_at_tested_concurrency=ideal))
        shared, scratch = [{row['file']: row['normalized_hits_sha256'] for row in mode['observations']} for mode in metrics]
        if shared != scratch:
            raise ValueError('Shared/scratch search evidence differs; do not adopt the storage change')
        if before != (database.stat().st_size, database.stat().st_mtime_ns):
            raise ValueError('Source database changed during the benchmark')
        write_json(args.out / 'benchmark.json', dict(status='PASS', hostname=os.uname().nodename,
            diamond_version=version, database=identity, database_bytes=before[0],
            jobs=args.jobs, threads=args.threads, staging_seconds=staging,
            cohort_chunks=selection['cohort_chunks'], modes=metrics,
            evidence_identical=True,
            estimate_limits='Sample extrapolation at tested concurrency; excludes download/build, scheduling and long-tail variation. Scratch is measured after shared reads and may benefit from warm caches.'))
    finally:
        shutil.rmtree(local)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('select')
    p.add_argument('--project', type=Path, required=True)
    p.add_argument('--trace', type=Path, required=True)
    p.add_argument('--samples', type=int, default=12)
    p.add_argument('--out', type=Path, required=True)
    b = commands.add_parser('run')
    b.add_argument('--queries', type=Path, required=True)
    b.add_argument('--bundle', type=Path, required=True)
    b.add_argument('--scratch', type=Path, required=True)
    b.add_argument('--jobs', type=int, default=8)
    b.add_argument('--threads', type=int, default=8)
    b.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'select':
        if args.samples < 1:
            parser.error('Samples must be positive')
        select(args.project.resolve(), args.trace.resolve(), args.samples, args.out)
    else:
        if min(args.jobs, args.threads) < 1:
            parser.error('Jobs and threads must be positive')
        def interrupted(signum, frame):
            raise KeyboardInterrupt('Benchmark interrupted')
        signal.signal(signal.SIGTERM, interrupted)
        benchmark(args)


if __name__ == '__main__':
    main()
