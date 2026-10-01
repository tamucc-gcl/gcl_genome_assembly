import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts/comparisons'))
from make_sv_benchmark import fixtures, reverse_complement


class BenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref, cls.cases, cls.truth = fixtures()

    def test_no_change_and_fragmentation_preserve_sequence(self):
        self.assertEqual(self.cases['control'], self.ref)
        for chrom, sequence in self.ref.items():
            self.assertEqual(''.join(self.cases['fragmented_control'][f'{chrom}_part{i}']
                                     for i in (1, 2, 3)), sequence)

    def test_truth_matches_inserted_or_inverted_sequence(self):
        for case, event, chrom, start, end, target, qstart, qend in self.truth:
            expected = self.ref[chrom][start:end]
            if event in ('INV', 'DUP_INVERTED_DISPERSED'):
                expected = reverse_complement(expected)
            self.assertEqual(self.cases[case][target][qstart:qend], expected)

    def test_translocation_conserves_total_bases_and_removes_source(self):
        query = self.cases['translocation']
        self.assertEqual(sum(map(len, query.values())), sum(map(len, self.ref.values())))
        self.assertEqual(query['chr1'], self.ref['chr1'][:100000] + self.ref['chr1'][120000:])


if __name__ == '__main__':
    unittest.main()
