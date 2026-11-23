import csv
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
from types import SimpleNamespace
import unittest

COMPARISONS = Path(__file__).resolve().parents[1] / 'scripts/comparisons'
sys.path.insert(0, str(COMPARISONS))
import contamination_staged_preflight as preflight


class PreflightTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.project = self.root / 'project'
        self.project.mkdir()
        self.trace = self.project / 'trace.tsv'
        self.out = self.project / 'preflight'
        rows = []
        for i in range(2):
            assembly = f'fixture_hap{i+1}'
            fasta_task = self.project / 'work' / f'{i:02d}' / 'abcdef012345'
            fasta_task.mkdir(parents=True)
            (fasta_task / f'{assembly}.decontaminated.fasta').write_text('>contig\n' + 'ACGT' * 300 + '\n')
            coverage_task = self.project / 'work' / f'{i+2:02d}' / 'abcdef012345'
            stats = coverage_task / f'{assembly}_mapping_stats'
            stats.mkdir(parents=True)
            (stats / 'coverage.txt').write_text('#rname\tstartpos\tendpos\tmeandepth\ncontig\t1\t1200\t20\n')
            rows.extend([
                dict(name=f'CONTIG_REFINEMENT:DECONTAMINATE_ASSEMBLY_CONTIG:FCS_CLEAN_GENOME ({assembly})',
                     hash=f'{i:02d}/abcdef', status='CACHED', exit='0'),
                dict(name=f'QC_PHASE:ASSEMBLY_QC_CONTIG_DECONTAM:MAPPING_QC ({assembly})',
                     hash=f'{i+2:02d}/abcdef', status='COMPLETED', exit='0')])
        with self.trace.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['name', 'hash', 'status', 'exit'], delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        self.args = SimpleNamespace(project=self.project, trace=self.trace, out=self.out, expected_assemblies=2,
                     workers=1, window_bp=10000, spacing_bp=1000000, short_contig_bp=100000,
                     gc_delta=.15, depth_fold=4)

    def tearDown(self):
        self.temporary.cleanup()

    def test_complete_census_reuses_identical_queries_and_excludes_fastas_from_archive(self):
        preflight.run(self.args)
        result = json.loads((self.out / 'preflight.json').read_text())
        self.assertEqual(result['totals']['protein_chunks'], 2)
        self.assertEqual(result['totals']['protein_query_bp'], 2400)
        self.assertEqual(result['alignment_jobs_launched'], 0)
        chunk = self.out / 'assemblies/fixture_hap1/queries/protein_00001.fasta'
        before = chunk.stat().st_mtime_ns
        preflight.run(self.args)
        self.assertEqual(chunk.stat().st_mtime_ns, before)
        self.assertTrue(all(row['reused'] == 'True' for row in preflight.table(self.out / 'summary.tsv')))
        reasons = preflight.table(self.out / 'selection_reasons.tsv')
        self.assertEqual(sum(int(row['query_bp']) for row in reasons), 2400)
        with tarfile.open(self.out.with_name('preflight.review.tar.gz')) as archive:
            self.assertTrue(all(not name.endswith('.fasta') for name in archive.getnames()))
            self.assertIn('summary.tsv', archive.getnames())

    def test_rejects_incomplete_cohort_mismatched_coverage_and_changed_queries(self):
        with self.assertRaisesRegex(ValueError, 'Expected 10'):
            preflight.discover(self.project, self.trace, 10)
        coverage = self.project / 'work/02/abcdef012345/fixture_hap1_mapping_stats/coverage.txt'
        original = coverage.read_text()
        coverage.write_text(original.replace('1200', '1199'))
        with self.assertRaisesRegex(ValueError, 'does not match'):
            preflight.run(self.args)
        coverage.write_text(original)
        preflight.run(self.args)
        chunk = self.out / 'assemblies/fixture_hap1/queries/protein_00001.fasta'
        chunk.write_text('>changed\nAAAA\n')
        with self.assertRaisesRegex(ValueError, 'Changed retained'):
            preflight.run(self.args)

    def test_two_process_cli_needs_no_alignment_executables(self):
        command = [sys.executable, str(COMPARISONS / 'contamination_staged_preflight.py'),
                   '--project', str(self.project), '--trace', str(self.trace), '--out', str(self.out),
                   '--expected-assemblies', '2', '--workers', '2']
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('no searches launched', result.stdout)
        self.assertEqual(json.loads((self.out / 'preflight.json').read_text())['status'], 'PASS')


if __name__ == '__main__':
    unittest.main()
