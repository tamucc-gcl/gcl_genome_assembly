"""Synthetic tests of exploratory bin-pair accounting; no sequencing files."""
import sys
import unittest
import tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/comparisons'))
from explore_chimera_region import local_contacts, assessment_fasta, sha


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
