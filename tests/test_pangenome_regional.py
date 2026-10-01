import csv
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
try:
    from pangenome_regional import aggregate, chromosome_label
except ModuleNotFoundError as exc:
    if exc.name != 'numpy':
        raise
    aggregate = None


@unittest.skipIf(aggregate is None, 'Regional adapter tests require numpy')
class RegionalTests(unittest.TestCase):
    def test_fragments_repeats_and_cross_chromosome_deduplication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'g.gfa').write_text(
                'S\t10\tAAAA\nS\t30\tCC\nS\t99\tGGG\n'
                'W\tA\t1\tchr1_1\t0\t6\t>10>30>10\n'
                'W\tA\t1\tchr1_1\t6\t8\t>30\n'
                'W\tA\t1\tchr2_1\t0\t3\t>99\n'
                'P\tA#2#chr1_1\t10+,99+\t*\nP\tB#1#chr1_1\t10+\t*\n')
            (root / 'groups').write_text('A#1#chr1_1\tA1\nA#1#chr2_1\tA1\nA#2#chr1_1\tA2\nB#1#chr1_1\tB1\n')
            (root / 'totals').write_text('panacus\ttable\nbp\ttotal\n10\t3\n30\t1\n99\t2\n')
            (root / 'hist').write_text('count\tbp\n0\t0\n1\t2\n2\t3\n3\t4\n')
            aggregate(root / 'g.gfa', root / 'totals', root / 'groups', root / 'hist', root / 'out')
            with (root / 'out/regional_coverage.tsv').open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            a1 = [r for r in rows if r['scope'] == 'haplotype' and r['haplotype'] == 'A1']
            self.assertEqual([int(r['bp']) for r in a1], [2, 3, 4])
            regional = [r for r in rows if r['scope'] == 'chromosome_haplotype' and r['haplotype'] == 'A1' and r['chromosome'] == 'chr1']
            self.assertEqual([int(r['bp']) for r in regional], [2, 0, 4])
            (root / 'totals').write_text('panacus\ttable\nbp\ttotal\n10\t1\n30\t1\n99\t2\n')
            with self.assertRaisesRegex(ValueError, 'reconcile'):
                aggregate(root / 'g.gfa', root / 'totals', root / 'groups', root / 'hist', root / 'bad')

    def test_names_do_not_split_composites(self):
        self.assertEqual(chromosome_label('chr7_2'), 'chr7')
        self.assertEqual(chromosome_label('chr7_1+chr12_1'), 'composite:chr7_1+chr12_1')
        self.assertEqual(chromosome_label('scaffold_7'), 'unplaced')


if __name__ == '__main__':
    unittest.main()
