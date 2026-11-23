import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch


HELPER = Path(__file__).resolve().parents[1] / 'py_scripts/contamination_workers.py'
if os.name == 'posix':
    spec = importlib.util.spec_from_file_location('contamination_workers', HELPER)
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)


@unittest.skipUnless(os.name == 'posix', 'Node workers use Linux file locks and CPU affinity')
class NodeWorkerTests(unittest.TestCase):
    def test_node_memory_records_sampled_available_ram(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'meminfo'
            path.write_text('MemTotal: 10000 kB\nMemAvailable: 8000 kB\n')
            monitor = worker.NodeMemory(path, interval=60)
            try:
                path.write_text('MemTotal: 10000 kB\nMemAvailable: 3000 kB\n')
                monitor.poll()
                path.write_text('MemTotal: 10000 kB\nMemAvailable: 7000 kB\n')
                result = monitor.finish()
                self.assertEqual(result['minimum_available_bytes'], 3000 * 1024)
                self.assertEqual(result['initial_available_bytes'], 8000 * 1024)
                self.assertEqual(result['last_available_bytes'], 7000 * 1024)
                self.assertEqual(result['total_bytes'], 10000 * 1024)
            finally:
                monitor.finish()

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        fake = self.bin / 'blastn'
        fake.write_text('''#!/usr/bin/env python3
import json, os, pathlib, sys, time
if '-version' in sys.argv:
    print('blastn: 2.16.0+'); raise SystemExit(0)
query=pathlib.Path(sys.argv[sys.argv.index('-query')+1]).read_text().splitlines()[0][1:]
with open(os.environ['CALLS'],'a') as out:
    out.write(query+'\\n'); out.flush()
if query=='B' and os.environ.get('SLOW_B'):
    pathlib.Path(os.environ['SLOW_B']).touch()
    time.sleep(1)
if query=='B' and os.environ.get('FAIL_ONCE'):
    marker=pathlib.Path(os.environ['FAIL_ONCE'])
    if not marker.exists():
        marker.touch(); raise SystemExit(1)
time.sleep(.05)
pathlib.Path(sys.argv[sys.argv.index('-out')+1]).write_text(query+'\\tmatch\\n')
''')
        fake.chmod(0o755)
        self.env = dict(os.environ, PATH=str(self.bin) + os.pathsep + os.environ['PATH'],
                        CALLS=str(self.root / 'calls.txt'))
        self.environment = patch.dict(os.environ, self.env)
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.addCleanup(self.temporary.cleanup)
        self.database = self.root / 'db' / 'nt'
        self.database.parent.mkdir()
        Path(str(self.database) + '.nal').write_text('DBLIST nt.000\n')
        (self.database.parent / 'nt.000.ndb').write_text('database fixture')
        self.specs = self.root / 'specs.json'
        specifications = []
        for assembly in ('A', 'B'):
            queries = self.root / assembly
            queries.mkdir()
            (queries / 'source.json').write_text(json.dumps(dict(assembly=assembly)))
            (queries / 'chunks.tsv').write_text('chunk\tqueries\tquery_bp\nchunk.fasta\t1\t4\n')
            (queries / 'chunk.fasta').write_text('>' + assembly + '\nACGT\n')
            specifications.append(dict(assembly=assembly, queries=str(queries)))
        self.specs.write_text(json.dumps(specifications))
        self.plan = self.root / 'plan.json'
        self.make_plan()

    def make_plan(self):
        worker.prepare(SimpleNamespace(database=self.database, specs=self.specs,
                                       state=self.root / 'state', out=self.plan, threads=1))
        return json.loads(self.plan.read_text())

    def start(self, identity=1, slots=1, **environment):
        return subprocess.Popen([sys.executable, str(HELPER), 'worker', '--plan', str(self.plan),
                                 '--worker', str(identity), '--slots', str(slots),
                                 '--scratch-root', str(self.root), '--shared-database',
                                 '--out', str(self.root / f'worker-{identity}.json')],
                                env=dict(self.env, **environment), stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL)

    def calls(self):
        path = self.root / 'calls.txt'
        return path.read_text().splitlines() if path.exists() else []

    def test_two_processes_claim_mixed_assemblies_once(self):
        processes = [self.start(1), self.start(2)]
        self.assertEqual([p.wait(timeout=15) for p in processes], [0, 0])
        self.assertCountEqual(self.calls(), ['A', 'B'])
        summaries = [json.loads((self.root / f'worker-{i}.json').read_text()) for i in (1, 2)]
        measured = [row['node_memory'] for row in summaries if row.get('searched')]
        self.assertTrue(measured)
        self.assertTrue(all(row['samples'] >= 2 and row['minimum_available_bytes'] > 0 for row in measured))
        # A fresh worker reuses the completed chunks without searching or staging again.
        self.assertEqual(self.start(3).wait(timeout=15), 0)
        self.assertCountEqual(self.calls(), ['A', 'B'])

    def test_failed_search_retains_and_reuses_successful_chunks(self):
        self.assertNotEqual(self.start(FAIL_ONCE=str(self.root / 'failed')).wait(timeout=15), 0)
        self.assertEqual(self.start(FAIL_ONCE=str(self.root / 'failed')).wait(timeout=15), 0)
        self.assertEqual(self.calls().count('A'), 1)
        self.assertEqual(self.calls().count('B'), 2)

    def test_interruption_reuses_committed_chunk(self):
        event = self.root / 'started'
        process = self.start(SLOW_B=str(event))
        deadline = time.monotonic() + 10
        while not event.exists() and time.monotonic() < deadline:
            time.sleep(.01)
        self.assertTrue(event.exists())
        process.send_signal(signal.SIGTERM)
        self.assertNotEqual(process.wait(timeout=10), 0)
        self.assertEqual(self.start().wait(timeout=15), 0)
        self.assertEqual(self.calls().count('A'), 1)
        self.assertEqual(self.calls().count('B'), 2)

    def test_corrupted_hits_are_recomputed(self):
        self.assertEqual(self.start().wait(timeout=15), 0)
        plan = json.loads(self.plan.read_text())
        store = worker.Results(plan)
        hits = store.directory(plan['tasks'][0]) / 'hits.tsv'
        hits.write_text('corruption')
        self.assertEqual(self.start().wait(timeout=15), 0)
        self.assertEqual(self.calls().count('A'), 2)
        self.assertEqual(self.calls().count('B'), 1)

    def test_hard_kill_releases_claim_and_keeps_committed_chunk(self):
        event = self.root / 'started'
        process = self.start(SLOW_B=str(event))
        deadline = time.monotonic() + 10
        while not event.exists() and time.monotonic() < deadline:
            time.sleep(.01)
        self.assertTrue(event.exists())
        process.kill()
        process.wait(timeout=10)
        # The separate fake BLAST process finishes without its killed parent committing a marker.
        time.sleep(1.2)
        self.assertEqual(self.start().wait(timeout=15), 0)
        self.assertEqual(self.calls().count('A'), 1)
        self.assertEqual(self.calls().count('B'), 2)

    def test_source_and_database_changes_invalidate_only_correct_keys(self):
        before = self.make_plan()
        (self.root / 'A/chunk.fasta').write_text('>A\nTGCA\n')
        changed = self.make_plan()
        self.assertNotEqual(before['tasks'][0]['key'], changed['tasks'][0]['key'])
        self.assertEqual(before['tasks'][1]['key'], changed['tasks'][1]['key'])
        (self.database.parent / 'nt.000.ndb').write_text('different database fixture')
        after = self.make_plan()
        self.assertNotEqual(changed['database_identity'], after['database_identity'])
        self.assertTrue(all(a['key'] != b['key'] for a, b in zip(changed['tasks'], after['tasks'])))

    def test_collector_releases_one_assembly_before_cohort_completion(self):
        collector = subprocess.Popen([sys.executable, str(HELPER), 'collect', '--plan', str(self.plan),
                                      '--assembly', 'A', '--out', str(self.root / 'collected'),
                                      '--poll-seconds', '.01', '--timeout', '10'])
        process = self.start(slots=2, SLOW_B=str(self.root / 'started'))
        self.assertEqual(collector.wait(timeout=15), 0)
        plan = json.loads(self.plan.read_text())
        self.assertIsNone(worker.Results(plan).completed(plan['tasks'][1]))
        self.assertEqual(process.wait(timeout=15), 0)
        self.assertEqual(len(list((self.root / 'collected').glob('*.hits.tsv'))), 1)

    def test_verified_database_copy(self):
        plan = self.make_plan()
        with patch.object(worker.shutil, 'disk_usage', return_value=SimpleNamespace(free=10**15, total=10**15)):
            local, seconds = worker.stage_database(plan, self.root / 'scratch-db')
        self.assertEqual(worker.digest(Path(str(local) + '.nal')),
                         worker.digest(Path(str(self.database) + '.nal')))
        self.assertEqual(worker.digest(local.parent / 'nt.000.ndb'),
                         worker.digest(self.database.parent / 'nt.000.ndb'))

    def test_changed_database_fails_worker_before_searching(self):
        (self.database.parent / 'nt.000.ndb').write_text('changed after plan')
        self.assertNotEqual(self.start().wait(timeout=15), 0)
        self.assertEqual(self.calls(), [])

    def test_space_limited_worker_releases_work_for_capable_worker(self):
        output = self.root / 'space-limited.json'
        args = SimpleNamespace(plan=self.plan, worker=1, slots=1, scratch_root=self.root,
                               shared_database=False, out=output)
        with patch.object(worker.shutil, 'disk_usage', return_value=SimpleNamespace(free=1, total=100)):
            worker.worker(args)
        skipped = json.loads(output.read_text())
        self.assertEqual(skipped['status'], 'skipped_insufficient_scratch')
        self.assertEqual(skipped['available_bytes'], 1)
        self.assertEqual(self.calls(), [])
        check_args = SimpleNamespace(plan=self.plan, summaries=[output], out=self.root / 'check.json')
        with self.assertRaisesRegex(RuntimeError, '2 chunks remain'):
            worker.check(check_args)
        self.assertEqual(self.start().wait(timeout=15), 0)
        worker.check(check_args)
        self.assertEqual(json.loads(check_args.out.read_text())['missing'], [])
        self.assertCountEqual(self.calls(), ['A', 'B'])


if __name__ == '__main__':
    unittest.main()
