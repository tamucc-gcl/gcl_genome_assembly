#!/usr/bin/env python3
"""Coordinate BLAST chunks across exclusive nodes with durable chunk completion records."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import tempfile
import threading
import time
import uuid


SEARCH = ['-task', 'megablast', '-evalue', '1e-10', '-culling_limit', '10',
          '-outfmt', '6 qseqid saccver staxids qstart qend bitscore evalue pident stitle']
SCHEMA = 'contamination-megablast-chunk-v1'


class InsufficientScratch(RuntimeError):
    def __init__(self, path, needed, usage):
        self.details = dict(hostname=os.uname().nodename, scratch=str(path),
                            required_bytes=needed, available_bytes=usage.free, total_bytes=usage.total)
        super().__init__(f'{self.details["hostname"]}: {path} requires {needed / 1024**3:.1f} GiB '
                         f'but has {usage.free / 1024**3:.1f} GiB available '
                         f'of {usage.total / 1024**3:.1f} GiB total')


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    with temporary.open('w') as handle:
        json.dump(value, handle, sort_keys=True, indent=2)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def blast_version():
    return subprocess.check_output(['blastn', '-version'], text=True).strip()


def database_files(prefix):
    """Keep BLAST volume/index/taxonomy files and their download provenance."""
    files = sorted(p for p in prefix.parent.iterdir() if p.is_file() and
                   (p.name.startswith(prefix.name + '.') or p.name.startswith('taxdb.') or
                    p.name.startswith('taxonomy4blast.') or p.name == 'database.json'))
    if not Path(str(prefix) + '.nal').is_file() and not Path(str(prefix) + '.ndb').is_file():
        raise ValueError(f'Missing BLAST nucleotide database: {prefix}')
    alias = Path(str(prefix) + '.nal')
    if alias.is_file():
        for line in alias.read_text().splitlines():
            if line.startswith('DBLIST '):
                if any('/' in token or '\\' in token for token in shlex.split(line)[1:]):
                    raise ValueError('Node-local staging requires alias volumes in the same database directory')
    return files


def database_identity(prefix):
    files = database_files(prefix)
    identity = [dict(name=p.name, size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns)
                for p in files]
    manifest = prefix.parent / 'database.json'
    value = dict(files=identity, manifest_sha256=digest(manifest) if manifest.is_file() else None)
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest(), identity


def prepare(args):
    prefix = args.database.resolve()
    identity, files = database_identity(prefix)
    version = blast_version()
    assemblies, by_assembly = [], []
    for specification in json.loads(args.specs.read_text()):
        assembly = specification['assembly']
        directory = Path(specification['queries']).resolve()
        source = json.loads((directory / 'source.json').read_text())
        if source['assembly'] != assembly:
            raise ValueError('Query manifest assembly does not match its workflow identity')
        if assembly in assemblies:
            raise ValueError('Duplicate assembly in search plan: ' + assembly)
        assemblies.append(assembly)
        import csv
        with (directory / 'chunks.tsv').open(newline='') as handle:
            chunks = list(csv.DictReader(handle, delimiter='\t'))
        tasks = []
        for chunk in chunks:
            query = (directory / chunk['chunk']).resolve(strict=True)
            if query.parent != directory:
                raise ValueError('Query chunk is outside its prepared directory')
            checksum = digest(query)
            key_data = dict(schema=SCHEMA, query_sha256=checksum, database=identity,
                            version=version, flags=SEARCH, threads=args.threads)
            key = hashlib.sha256(json.dumps(key_data, sort_keys=True).encode()).hexdigest()
            tasks.append(dict(assembly=assembly, query=str(query), name=chunk['chunk'],
                              query_sha256=checksum, key=key))
        by_assembly.append(tasks)
    # Round-robin assembly chunks so the first assembly cannot occupy the whole queue.
    tasks = [chunks[i] for i in range(max(map(len, by_assembly), default=0))
             for chunks in by_assembly if i < len(chunks)]
    atomic_json(args.out, dict(schema=SCHEMA, database=str(prefix), database_identity=identity,
                             database_files=files, version=version, flags=SEARCH, threads=args.threads,
                             state=str(args.state.resolve()), assemblies=assemblies, tasks=tasks))


class Results:
    def __init__(self, plan):
        self.plan = plan
        self.root = Path(plan['state'])
        self.validated = {}

    def directory(self, task):
        return self.root / task['key'][:2] / task['key']

    def completed(self, task):
        directory = self.directory(task)
        marker, hits = directory / 'completion.json', directory / 'hits.tsv'
        try:
            record = json.loads(marker.read_text())
            stat = hits.stat()
            signature = (stat.st_size, stat.st_mtime_ns, record.get('hits_sha256'))
            if record.get('key') != task['key'] or record.get('query_sha256') != task['query_sha256']:
                return None
            if record.get('hits_bytes') != stat.st_size:
                return None
            if self.validated.get(task['key']) != signature:
                if digest(hits) != record.get('hits_sha256'):
                    return None
                self.validated[task['key']] = signature
            return hits
        except (OSError, ValueError):
            return None


def stage_database(plan, directory, stopped=None):
    prefix = Path(plan['database'])
    current, _ = database_identity(prefix)
    if current != plan['database_identity']:
        raise ValueError('The source nucleotide database changed after search planning')
    needed = sum(row['size'] for row in plan['database_files'])
    usage = shutil.disk_usage(directory.parent)
    if usage.free < needed + 10 * 1024 ** 3:
        raise InsufficientScratch(directory.parent, needed + 10 * 1024 ** 3, usage)
    directory.mkdir()
    start = time.monotonic()
    for row in plan['database_files']:
        source, destination = prefix.parent / row['name'], directory / row['name']
        checksum = hashlib.sha256()
        with source.open('rb') as reader, destination.open('wb') as writer:
            for block in iter(lambda: reader.read(8 * 1024 ** 2), b''):
                if stopped is not None and stopped.is_set():
                    raise RuntimeError('Worker interrupted during database staging')
                checksum.update(block)
                writer.write(block)
        if destination.stat().st_size != row['size'] or digest(destination) != checksum.hexdigest():
            raise ValueError('Scratch database copy failed verification: ' + row['name'])
    if database_identity(prefix)[0] != plan['database_identity']:
        raise ValueError('The source nucleotide database changed during staging')
    return directory / prefix.name, time.monotonic() - start


def worker(args):
    plan = json.loads(args.plan.read_text())
    if blast_version() != plan['version']:
        raise ValueError('Worker BLAST version differs from its search plan')
    store = Results(plan)
    pending = [task for task in plan['tasks'] if not store.completed(task)]
    if not pending:
        atomic_json(args.out, dict(worker=args.worker, searched=0, reused=len(plan['tasks'])))
        return
    cpus = sorted(os.sched_getaffinity(0))
    if len(cpus) < args.slots * plan['threads']:
        raise ValueError('Worker allocation has fewer CPUs than slots multiplied by threads')
    local = Path(tempfile.mkdtemp(prefix='contamination-nt-', dir=args.scratch_root))
    stopped = threading.Event()
    children = {}
    outcomes = []

    def terminate(signum, frame):
        stopped.set()
        for process in tuple(children.values()):
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    old_handlers = {s: signal.signal(s, terminate) for s in (signal.SIGTERM, signal.SIGINT)}
    try:
        if args.shared_database:
            if database_identity(Path(plan['database']))[0] != plan['database_identity']:
                raise ValueError('The source nucleotide database changed after planning')
            prefix, staging_seconds = Path(plan['database']), 0
        else:
            try:
                prefix, staging_seconds = stage_database(plan, local / 'database', stopped)
            except InsufficientScratch as error:
                print(f'Worker {args.worker} skipped: {error}', flush=True)
                atomic_json(args.out, dict(worker=args.worker, searched=0,
                                           status='skipped_insufficient_scratch', **error.details))
                return
        print(f'Worker {args.worker} ready on {os.uname().nodename}; staging {staging_seconds:.1f}s', flush=True)

        def slot(index):
            cursor = ((args.worker - 1) * args.slots + index) % len(pending)
            assigned = cpus[index * plan['threads']:(index + 1) * plan['threads']]
            while not stopped.is_set():
                claimed = False
                for offset in range(len(pending)):
                    position = (cursor + offset) % len(pending)
                    task = pending[position]
                    if store.completed(task):
                        continue
                    directory = store.directory(task)
                    directory.mkdir(parents=True, exist_ok=True)
                    with (directory / 'claim.lock').open('a') as lock:
                        try:
                            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        except BlockingIOError:
                            continue
                        if store.completed(task):
                            continue
                        claimed = True
                        cursor = (position + 1) % len(pending)
                        attempt = uuid.uuid4().hex
                        query = local / f'{attempt}.fasta'
                        hits = local / f'{attempt}.hits.tsv'
                        shutil.copyfile(task['query'], query)
                        if digest(query) != task['query_sha256']:
                            raise ValueError('Query sequence changed after planning: ' + task['query'])
                        start = time.monotonic()
                        command = ['taskset', '-c', ','.join(map(str, assigned)), 'blastn',
                                   '-db', str(prefix), '-query', str(query), '-out', str(hits),
                                   '-num_threads', str(plan['threads']), *plan['flags']]
                        stderr = directory / f'attempt-{attempt}.stderr.txt'
                        with stderr.open('w') as error:
                            process = subprocess.Popen(command, stderr=error, start_new_session=True,
                                                       env=dict(os.environ, BLASTDB=str(prefix.parent)))
                            children[index] = process
                            if stopped.is_set():
                                terminate(None, None)
                            status = process.wait()
                            children.pop(index, None)
                        if status:
                            raise RuntimeError(f'BLAST failed with exit {status}; see {stderr}')
                        if stopped.is_set():
                            break
                        temporary = directory / f'hits-{attempt}.tmp'
                        shutil.copyfile(hits, temporary)
                        checksum = digest(temporary)
                        os.replace(temporary, directory / 'hits.tsv')
                        record = dict(key=task['key'], query_sha256=task['query_sha256'],
                                      database_identity=plan['database_identity'], version=plan['version'],
                                      flags=plan['flags'], threads=plan['threads'],
                                      hits_sha256=checksum, hits_bytes=hits.stat().st_size,
                                      seconds=time.monotonic() - start, worker=args.worker,
                                      hostname=os.uname().nodename)
                        atomic_json(directory / 'completion.json', record)
                        outcomes.append(record)
                        query.unlink()
                        hits.unlink()
                        print(f'Completed {task["assembly"]}/{task["name"]}', flush=True)
                    break
                if not claimed:
                    # Every unfinished chunk is already owned by another active slot/node.
                    # Release idle node resources instead of waiting for the cohort's last chunk.
                    return

        def guarded_slot(index):
            try:
                slot(index)
            except BaseException:
                terminate(None, None)
                raise

        with ThreadPoolExecutor(max_workers=args.slots) as pool:
            list(pool.map(guarded_slot, range(args.slots)))
        if stopped.is_set():
            raise RuntimeError('Worker interrupted; completed chunks have been retained')
        atomic_json(args.out, dict(worker=args.worker, searched=len(outcomes),
                                   staging_seconds=staging_seconds, completed=outcomes))
    finally:
        for signum, handler in old_handlers.items():
            signal.signal(signum, handler)
        shutil.rmtree(local)


def collect(args):
    plan = json.loads(args.plan.read_text())
    if args.assembly not in plan['assemblies']:
        raise ValueError('Assembly missing from node search plan')
    tasks = [task for task in plan['tasks'] if task['assembly'] == args.assembly]
    store = Results(plan)
    args.out.mkdir()
    remaining = list(enumerate(tasks))
    deadline = time.monotonic() + args.timeout
    while remaining:
        missing = []
        for index, task in remaining:
            hits = store.completed(task)
            if hits:
                target = args.out / f'{index:06d}_{task["name"]}.hits.tsv'
                shutil.copyfile(hits, target)
                if digest(target) != digest(hits):
                    raise ValueError('Collected output differs from completed search')
            else:
                missing.append((index, task))
        remaining = missing
        if remaining:
            if time.monotonic() >= deadline:
                raise TimeoutError('Timed out waiting for assembly chunks: ' + args.assembly)
            time.sleep(args.poll_seconds)


def check(args):
    plan = json.loads(args.plan.read_text())
    store = Results(plan)
    missing = [dict(assembly=task['assembly'], name=task['name'])
               for task in plan['tasks'] if not store.completed(task)]
    summaries = [json.loads(path.read_text()) for path in args.summaries]
    skipped = [row for row in summaries if row.get('status') == 'skipped_insufficient_scratch']
    atomic_json(args.out, dict(chunks=len(plan['tasks']), missing=missing, skipped_workers=skipped))
    if missing:
        raise RuntimeError(f'{len(missing)} chunks remain after all node workers exited. '
                           f'{len(skipped)} workers lacked scratch space. See {args.out}.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('plan')
    for name in ('specs', 'database', 'state', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--threads', type=int, default=8)
    w = commands.add_parser('worker')
    w.add_argument('--plan', type=Path, required=True)
    w.add_argument('--worker', type=int, required=True)
    w.add_argument('--slots', type=int, default=8)
    w.add_argument('--scratch-root', type=Path, default=Path('/scratch'))
    w.add_argument('--shared-database', action='store_true', help='Controlled testing without staging')
    w.add_argument('--out', type=Path, required=True)
    c = commands.add_parser('collect')
    c.add_argument('--plan', type=Path, required=True)
    c.add_argument('--assembly', required=True)
    c.add_argument('--out', type=Path, required=True)
    c.add_argument('--poll-seconds', type=float, default=10)
    c.add_argument('--timeout', type=float, default=4 * 86400)
    v = commands.add_parser('check')
    v.add_argument('--plan', type=Path, required=True)
    v.add_argument('--summaries', type=Path, nargs='+', required=True)
    v.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan' and args.threads < 1:
        parser.error('Threads must be positive')
    if args.command == 'worker' and (args.slots < 1 or args.worker < 1):
        parser.error('Slots and worker identity must be positive')
    if args.command == 'collect' and (args.poll_seconds <= 0 or args.timeout <= 0):
        parser.error('Collector polling interval and timeout must be positive')
    dict(plan=prepare, worker=worker, collect=collect, check=check)[args.command](args)


if __name__ == '__main__':
    main()
