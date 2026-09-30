"""Synthetic tests of exploratory bin-pair accounting; no sequencing files."""
import sys
import unittest
import tempfile
import gzip
import json
from unittest.mock import patch
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/comparisons'))
from explore_chimera_region import local_contacts, assessment_fasta, sha, retained_file, main


class StagedFastaTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_flat_renamed_input_without_extension(self):
        path = self.root/'assessment'
        path.write_text('>s\nACGT\n')
        self.assertEqual(assessment_fasta([self.root], sha(path)), path)

    def test_nested_staging_with_index(self):
        nested = self.root/'assessment'
        nested.mkdir()
        path = nested/'assembly.fasta'
        path.write_text('>s\nACGT\n')
        (nested/'assembly.fasta.fai').write_text('s\t4\t3\t4\t5\n')
        self.assertEqual(assessment_fasta([self.root], sha(path)), path)

    def test_fallback_to_evidence_task_rejects_wrong_stage(self):
        joins, evidence = self.root/'joins', self.root/'evidence'
        joins.mkdir()
        evidence.mkdir()
        (joins/'wrong.fa').write_text('>s\nAAAA\n')
        path = evidence/'current.fa'
        path.write_text('>s\nACGT\n')
        self.assertEqual(assessment_fasta([joins, evidence], sha(path)), path)

    def test_missing_match_fails_with_inventory(self):
        (self.root/'wrong.fa').write_text('>s\nAAAA\n')
        with self.assertRaisesRegex(ValueError, 'Task entries:.*wrong.fa'):
            assessment_fasta([self.root], '0'*64)

    def test_identical_copies_are_not_ambiguous(self):
        first = self.root/'a.fa'
        first.write_text('>s\nACGT\n')
        (self.root/'b.fa').write_text(first.read_text())
        self.assertEqual(assessment_fasta([self.root], sha(first)), first)

    def test_explicit_wrong_stage_is_rejected(self):
        path = self.root/'wrong.fa'
        path.write_text('>s\nAAAA\n')
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            assessment_fasta([], '0'*64, explicit=path)

    def test_scratch_input_falls_back_to_published_output(self):
        path = self.root/'published.tsv'
        path.write_text('retained output\n')
        self.assertEqual(retained_file([], path, 'test input'), path)

    def test_missing_published_input_fails_explicitly(self):
        with self.assertRaisesRegex(ValueError, 'Missing test input'):
            retained_file([], self.root/'missing.tsv', 'test input')

    def test_end_to_end_with_only_declared_task_outputs(self):
        joins, evidence, results = (self.root/name for name in ('joins', 'evidence', 'results'))
        for path in (joins, evidence, results):
            path.mkdir()
        # No staged inputs and no projected pairs survive in either task directory.
        fasta = results/'pre_finish.fa'
        fasta.write_text('>s\n'+'A'*5000000+'\n')
        calls = ('assembly\tscaffold\tassembly_sha256\tjoin_scope\n'
                 'a\ts\t'+sha(fasta)+'\tround2_only\n')
        (joins/'a.chimeric_joins.tsv').write_text(calls)
        review = results/'assembly/chimeras/older_joins'
        review.mkdir(parents=True)
        (review/'a.review_joins.tsv').write_text(calls)
        agp = results/'last.agp'
        agp.write_text('s\t1\t5000000\t1\tW\tcomponent\t1\t5000000\t+\n')
        source = results/'source.pairs.gz'
        with gzip.open(source, 'wt') as f:
            f.write('r1\tcomponent\t2400001\tcomponent\t2600001\t+\t-\tUU\n')
        alignments = results/'assembly/harmonization'
        alignments.mkdir(parents=True)
        with gzip.open(alignments/'a.ref.paf.gz', 'wt') as f:
            f.write('s\t5000000\t0\t5000000\t+\tref\t5000000\t0\t5000000\t5000000\t5000000\t60\n')
        log = self.root/'run.log'
        log.write_text('')
        out = self.root/'diagnostic'
        args = ['explore_chimera_region.py', '--joins-work', str(joins),
                '--evidence-work', str(evidence), '--results', str(results),
                '--work-root', str(self.root), '--log', str(log), '--assembly', 'a',
                '--scaffold', 's', '--start', '2000000', '--end', '3000000',
                '--out', str(out), '--assessment-fasta', str(fasta),
                '--last-agp', str(agp), '--source-pairs', str(source)]
        with patch.object(sys, 'argv', args):
            main()
        self.assertTrue(json.loads((out/'provenance.json').read_text())['pairs_rebuilt'])
        self.assertTrue((out/'local_hic.tsv').is_file())
        self.assertIn('agp_fasta_check\tpassed', (out/'coordinate_audit.tsv').read_text())
        self.assertIn('r1\ts\t2400001\ts\t2600001',
                      (out/'rebuilt_pairs/a.s.pairs').read_text())


class RegionalContacts(unittest.TestCase):
    def test_uniform_contacts_have_equal_means(self):
        matrix = {(x,y): 3 for x in range(8) for y in range(x+1,8)}
        totals, opportunities = local_contacts(matrix, 4, 4, 1, 4)
        for label in ('cross', 'left', 'right'):
            self.assertEqual(totals[label]/opportunities[label], 3)

    def test_missing_cross_contacts_do_not_remove_opportunities(self):
        matrix = {(x,y): 3 for x in range(8) for y in range(x+1,8) if (x < 4) == (y < 4)}
        totals, opportunities = local_contacts(matrix, 4, 4, 1, 4)
        self.assertEqual(totals['cross'], 0)
        self.assertGreater(opportunities['cross'], 0)
        self.assertEqual(totals['left']/opportunities['left'], 3)

    def test_diagonal_and_outside_distance_are_excluded(self):
        totals, _ = local_contacts({(4,4): 100, (0,7): 100, (3,4): 2}, 4, 4, 1, 2)
        self.assertEqual(sum(totals.values()), 2)


if __name__ == '__main__':
    unittest.main()
