import json
from pathlib import Path
import random
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from contamination_review import prepare, records, write_table
from contamination_evidence import table, select, summarize, rebuild_report
from contamination_staged import sample, route, coding_score


class StagedTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.fasta = self.root / 'source.fasta'
        rng = random.Random(26)
        self.sequences = dict(host=''.join(rng.choices('ACGT', k=200000)),
                              short=''.join(rng.choices('ACGT', k=3000)),
                              depth=''.join(rng.choices('ACGT', k=10000)), ambiguous='N' * 2000)
        self.fasta.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in self.sequences.items()))
        self.metrics = self.root / 'metrics'
        self.metrics.mkdir()
        prepare(SimpleNamespace(fasta=self.fasta, assembly='fixture', stage='post_fcs'), self.metrics)
        self.coverage = self.root / 'coverage.tsv'
        write_table(self.coverage, [dict(zip(['#rname', 'startpos', 'endpos', 'meandepth'],
                            [name, 1, len(seq), 300 if name == 'depth' else 30]))
                            for name, seq in self.sequences.items()], ['#rname', 'startpos', 'endpos', 'meandepth'])
        self.queries = self.root / 'queries'
        self.args = SimpleNamespace(fasta=self.fasta, metrics=self.metrics, coverage=self.coverage,
             out=self.queries, window_bp=10000, spacing_bp=50000, chunk_bp=50000,
             short_contig_bp=5000, gc_delta=.15, depth_fold=4)
        self.taxonomy = self.root / 'taxonomy'
        self.taxonomy.mkdir()
        (self.taxonomy / 'nodes.dmp').write_text('1 | 1 | no rank |\n2 | 1 | phylum |\n3 | 1 | phylum |\n10 | 2 | species |\n20 | 3 | species |\n')
        (self.taxonomy / 'names.dmp').write_text('2 | Chordata | | scientific name |\n3 | Arthropoda | | scientific name |\n')

    def tearDown(self):
        self.tmp.cleanup()

    def routed(self):
        sample(self.args)
        rows = table(self.queries / 'windows.tsv')
        self.primary = [row for row in rows if row['status'] == 'search']
        self.hits = self.root / 'protein.hits.tsv'
        self.hits.write_text(''.join(f"{row['query']}\tsubject\t{20 if row['sequence'] == 'short' else 10}\t1\t1500\t100\t1e-30\t90\ttitle\n"
                    for row in self.primary if row['sequence'] != 'depth'))
        self.nt = self.root / 'nt'
        route(SimpleNamespace(fasta=self.fasta, prepared=self.queries, hits=[self.hits],
              taxdump=self.taxonomy, taxid=10, host_rank='phylum', out=self.nt))
        return rows

    def test_sampling_and_routing_preserve_coordinates_and_bound_work(self):
        rows = self.routed()
        host_samples = [row for row in self.primary if row['sequence'] == 'host']
        self.assertEqual(len(host_samples), 4)
        self.assertEqual({int(row['start0']) // 50000 for row in host_samples}, {0, 1, 2, 3})
        self.assertEqual(len(self.primary), 6)
        self.assertEqual({row['sequence'] for row in table(self.nt / 'windows.tsv')}, {'short', 'depth'})
        self.assertTrue(any(row['status'] == 'unsearched' for row in rows))
        self.assertEqual(sum(int(row['end0']) - int(row['start0']) for row in rows), 215000)
        for directory in (self.queries, self.nt):
            lookup = {row['query']: row for row in table(directory / 'windows.tsv')}
            for chunk in directory.glob('protein_*.fasta'):
                found = dict(records(chunk))
                self.assertLessEqual(sum(map(len, found.values())), 50000)
                for query, seq in found.items():
                    row = lookup[query]
                    self.assertEqual(seq, self.sequences[row['sequence']][int(row['start0']):int(row['end0'])])

    def test_flanks_are_searched_by_nt_without_claiming_protein_search(self):
        self.routed()
        # No protein hits: all representative windows and their adjacent windows route to nt.
        self.hits.write_text('')
        self.nt = self.root / 'nt-nohit'
        route(SimpleNamespace(fasta=self.fasta, prepared=self.queries, hits=[],
              taxdump=self.taxonomy, taxid=10, host_rank='phylum', out=self.nt))
        selected = self.root / 'selected'
        select(SimpleNamespace(fasta=self.fasta, prepared=self.queries, hits=[], megablast_hits=[],
               megablast=True, megablast_plan=self.nt / 'windows.tsv', followup_max_bp=100000,
               taxdump=self.taxonomy, taxid=10, host_rank='phylum', out=selected))
        windows = table(selected / 'windows.tsv')
        self.assertTrue(any(row['protein_evidence'] == 'unsearched' and row['megablast_evidence'] == 'no-hit'
                            for row in windows))
        self.assertTrue(any(row['megablast_evidence'] == 'unsearched' for row in windows))
        self.assertTrue(any(row['protein_evidence'] == 'insufficient-acgt' for row in windows))

    def test_report_and_rebuild_disclose_sampling_without_removal(self):
        self.routed()
        selected = self.root / 'selected'
        select(SimpleNamespace(fasta=self.fasta, prepared=self.queries, hits=[self.hits], megablast_hits=[],
               megablast=True, megablast_plan=self.nt / 'windows.tsv', followup_max_bp=100000,
               taxdump=self.taxonomy, taxid=10, host_rank='phylum', out=selected))
        database = self.root / 'protein_database.json'
        database.write_text('{}')
        review = self.root / 'review'
        summarize(SimpleNamespace(prepared=self.queries, selected=selected, metrics=self.metrics,
                  database=database, coverage=self.coverage, coverage_source='fixture',
                  nucleotide_hits=[], targeted=False, nucleotide_db='nt', nucleotide_manifest=None,
                  taxdump=self.taxonomy, taxid=10, host_rank='phylum', out=review))
        metadata = json.loads((review / 'evidence.json').read_text())
        self.assertEqual(metadata['search_mode'], 'staged')
        self.assertEqual(metadata['protein_searched_bp'], 53000)
        self.assertEqual(metadata['megablast_searched_bp'], 13000)
        self.assertEqual(metadata['sampling']['depth_baseline'], 30)
        self.assertIn('sampling_windows.tsv', (review / 'README.md').read_text())
        self.assertTrue(all(row['selected'] == 'NO' for row in table(review / 'review.tsv')))
        rebuilt = self.root / 'rebuilt'
        rebuild_report(SimpleNamespace(previous=review, fasta=self.fasta, out=rebuilt))
        self.assertEqual(json.loads((rebuilt / 'evidence.json').read_text())['sampling'], metadata['sampling'])
        self.assertEqual((review / 'sampling_windows.tsv').read_bytes(), (rebuilt / 'sampling_windows.tsv').read_bytes())

    def test_all_ambiguous_emits_no_search_chunks(self):
        self.fasta.write_text('>ambiguous\n' + 'N' * 2000 + '\n')
        # Recreate source-bound metrics for this input.
        self.metrics = self.root / 'ambiguous-metrics'
        self.metrics.mkdir()
        prepare(SimpleNamespace(fasta=self.fasta, assembly='fixture', stage='post_fcs'), self.metrics)
        self.args.metrics, self.args.coverage = self.metrics, None
        sample(self.args)
        self.assertEqual(table(self.queries / 'chunks.tsv'), [])
        self.assertEqual(json.loads((self.queries / 'source.json').read_text())['search_query_bp'], 0)

    def test_gc_outlier_and_missing_coverage_are_explicit(self):
        self.fasta.write_text('>host\n' + 'ATGC' * 50000 + '\n>gc\n' + 'G' * 10000 + '\n')
        metrics = self.root / 'gc-metrics'
        metrics.mkdir()
        prepare(SimpleNamespace(fasta=self.fasta, assembly='fixture', stage='post_fcs'), metrics)
        self.args.metrics, self.args.coverage = metrics, None
        sample(self.args)
        gc = [row for row in table(self.queries / 'windows.tsv') if row['sequence'] == 'gc']
        self.assertIn('gc-outlier', gc[0]['selection_reason'])
        self.assertFalse(json.loads((self.queries / 'source.json').read_text())['sampling']['coverage_available'])

    def test_coding_score_does_not_promote_ambiguous_runs(self):
        self.assertGreater(coding_score('ATG' * 200), coding_score('N' * 10000))

    def test_source_identity_and_plan_coordinates_fail_closed(self):
        self.routed()
        plan = table(self.nt / 'windows.tsv')
        plan[0]['start0'] = int(plan[0]['start0']) + 1
        write_table(self.nt / 'windows.tsv', plan, list(plan[0]))
        with self.assertRaisesRegex(ValueError, 'coordinates differ'):
            select(SimpleNamespace(fasta=self.fasta, prepared=self.queries, hits=[self.hits], megablast_hits=[],
                   megablast=True, megablast_plan=self.nt / 'windows.tsv', followup_max_bp=100000,
                   taxdump=self.taxonomy, taxid=10, host_rank='phylum', out=self.root / 'invalid'))
        self.fasta.write_text('>different\n' + 'A' * 1000 + '\n')
        with self.assertRaisesRegex(ValueError, 'differs from metrics source'):
            sample(SimpleNamespace(**dict(vars(self.args), out=self.root / 'changed')))


if __name__ == '__main__':
    unittest.main()
