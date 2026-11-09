"""Known-sequence tests for manual contamination removal, independent of a cohort."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1]/'py_scripts'/'contamination_review.py'
spec = importlib.util.spec_from_file_location('contamination_review', SCRIPT)
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


class ContaminationReviewTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root/'source.fa'
        self.source.write_text('>a original description\nACGTacgtNN\n>b\nTTTTCCCC\n')
        self.checksum = review.digest(self.source)

    def row(self, name='a', lo=2, hi=5, identity='manual1', action='REMOVE_INTERVAL'):
        return dict(zip(review.FIELDS, [identity, 'sample_hap1', 'post_fcs', self.checksum,
                    name, action, str(lo), str(hi), 'YES', 'reviewer', 'foreign sequence']))

    def run_tool(self, command, rows=None, success=True):
        dest = self.root/('result_'+str(len(list(self.root.iterdir()))))
        cmd = [sys.executable, str(SCRIPT), command, '--fasta', str(self.source),
               '--assembly', 'sample_hap1', '--stage', 'post_fcs', '--out', str(dest)]
        if command == 'apply':
            sheet = self.root/'actions.tsv'
            review.write_table(sheet, rows, list(rows[0]) if rows else review.FIELDS)
            cmd += ['--review', str(sheet)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(dest.exists())
            self.assertFalse(list(self.root.glob('.contamination-*')))
        return dest

    def test_internal_removal_splits_and_reconstructs_exact_bases(self):
        dest = self.run_tool('apply', [self.row()])
        retained = dict(review.records(dest/'retained.fasta'))
        removed = dict(review.records(dest/'removed.fasta'))
        self.assertEqual(retained, {'a__retained_0_2':'AC', 'a__retained_5_10':'cgtNN', 'b':'TTTTCCCC'})
        self.assertEqual(removed, {'a__removed_2_5':'GTa'})
        with open(dest/'coordinate_map.tsv') as handle:
            lift = list(csv.DictReader(handle, delimiter='\t'))
        pieces = {**retained, **removed}
        reconstructed = ''.join(pieces[row['output_sequence']] for row in lift if row['sequence']=='a')
        self.assertEqual(reconstructed, 'ACGTacgtNN')
        data = json.loads((dest/'verification.json').read_text())
        self.assertEqual((data['source_bp'], data['retained_bp'], data['removed_bp']), (18, 15, 3))

    def test_terminal_intervals_trim_without_discarding_other_sequence(self):
        dest = self.run_tool('apply', [self.row(lo=0, hi=2), self.row(lo=8, hi=10, identity='manual2')])
        self.assertEqual(dict(review.records(dest/'retained.fasta'))['a__retained_2_8'], 'GTacgt')

    def test_whole_sequence_removal(self):
        dest = self.run_tool('apply', [self.row(lo=0, hi=10, action='REMOVE_SEQUENCE')])
        self.assertEqual(dict(review.records(dest/'retained.fasta')), {'b':'TTTTCCCC'})

    def test_no_selection_preserves_fasta_bytes(self):
        row = self.row(); row['selected']='NO'
        dest = self.run_tool('apply', [row])
        self.assertEqual((dest/'retained.fasta').read_bytes(), self.source.read_bytes())
        self.assertEqual((dest/'removed.fasta').stat().st_size, 0)

    def test_empty_sheet_preserves_fasta_bytes(self):
        dest = self.run_tool('apply', [])
        self.assertEqual((dest/'retained.fasta').read_bytes(), self.source.read_bytes())

    def test_overlapping_and_duplicate_removals_rejected(self):
        self.run_tool('apply', [self.row(), self.row(lo=4, hi=8, identity='manual2')], success=False)
        self.run_tool('apply', [self.row(), self.row(identity='manual2')], success=False)

    def test_adjacent_removals_are_valid(self):
        dest = self.run_tool('apply', [self.row(lo=2, hi=5), self.row(lo=5, hi=7, identity='manual2')])
        retained = dict(review.records(dest/'retained.fasta'))
        self.assertEqual(retained['a__retained_7_10'], 'tNN')
        self.assertEqual(len(retained), 3)

    def test_checksum_stage_and_review_required(self):
        for key, value in [('source_sha256','bad'), ('stage','final'), ('reviewer',''),
                           ('reason',''), ('selected','yes'), ('action','CUT')]:
            with self.subTest(key=key):
                row = self.row(); row[key]=value
                self.run_tool('apply', [row], success=False)

    def test_bounds_and_whole_sequence_semantics(self):
        for lo, hi, action in [(-1, 5,'REMOVE_INTERVAL'), (2,11,'REMOVE_INTERVAL'),
                              (2,2,'REMOVE_INTERVAL'), (2,5,'REMOVE_SEQUENCE')]:
            self.run_tool('apply', [self.row(lo=lo, hi=hi, action=action)], success=False)

    def test_no_entire_assembly_deletion(self):
        self.run_tool('apply', [self.row(lo=0,hi=10), self.row(name='b',lo=0,hi=8,identity='manual2')], success=False)

    def test_prepare_metrics_and_source_bound_unselected_template(self):
        dest = self.run_tool('prepare')
        rows = review.read_actions(dest/'review.tsv')
        self.assertEqual([row['selected'] for row in rows], ['NO','NO'])
        self.assertTrue(all(row['source_sha256']==self.checksum for row in rows))
        with open(dest/'sequence_metrics.tsv') as handle:
            metrics = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(metrics[0]['n_bp'], '2')
        self.assertEqual(float(metrics[0]['gc_fraction']), 0.5)

    def test_output_identifier_collision_fails_before_publishing(self):
        self.source.write_text('>a\nACGTacgtNN\n>a__retained_0_2\nTTTT\n')
        self.checksum = review.digest(self.source)
        self.run_tool('apply', [self.row()], success=False)

    def test_duplicate_fasta_identifiers_and_invalid_bases(self):
        for content in ['>a\nACGT\n>a\nTTTT\n', '>a\nACGT-\n', '>a\n']:
            self.source.write_text(content)
            self.run_tool('prepare', success=False)

    def test_review_table_can_include_evidence_columns(self):
        row = self.row(); row['evidence']='protein and nucleotide'
        dest = self.run_tool('apply', [row])
        self.assertIn('evidence', (dest/'applied_actions.tsv').read_text().splitlines()[0])

    def test_search_windows_keep_coordinates_and_bound_chunk_size(self):
        from argparse import Namespace
        output = self.root/'windows'; output.mkdir()
        args = Namespace(fasta=self.source, assembly='sample_hap1', stage='post_fcs',
                         window_bp=4, chunk_bp=8, min_acgt_bp=2)
        review.search_input(args, output)
        with open(output/'windows.tsv') as handle:
            windows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(windows), 5)
        self.assertEqual(windows[2]['status'], 'insufficient_acgt')
        self.assertEqual((windows[2]['start0'], windows[2]['end0']), ('8','10'))
        data = json.loads((output/'source.json').read_text())
        self.assertEqual((data['source_bp'],data['search_query_bp'],data['chunks']), (18,16,2))
        queries = {}
        for chunk in output.glob('protein_*.fasta'):
            content = dict(review.records(chunk))
            self.assertLessEqual(sum(map(len, content.values())), 8)
            queries.update(content)
        originals = {header.split()[0]:seq for header,seq in review.records(self.source)}
        for row in windows:
            if row['status']=='search':
                self.assertEqual(queries[row['query']], originals[row['sequence']][int(row['start0']):int(row['end0'])])

    def test_all_ambiguous_search_input_explicitly_has_no_tasks(self):
        from argparse import Namespace
        self.source.write_text('>ambiguous\nNNNNNNNN\n')
        output = self.root/'windows'; output.mkdir()
        review.search_input(Namespace(fasta=self.source, assembly='sample',stage='post_fcs',
            window_bp=4,chunk_bp=8,min_acgt_bp=1), output)
        self.assertFalse(list(output.glob('protein_*.fasta')))
        self.assertEqual(json.loads((output/'source.json').read_text())['search_queries'],0)


if __name__ == '__main__':
    unittest.main()
