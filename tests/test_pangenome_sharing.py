import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from pangenome_sharing import prepare, read_hist, summarize, validate_growth


class SharingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.ledger = self.root / 'identity.tsv'
        self.ledger.write_text('id\tsample\tgraph_name\nA1\tA\tREF_A1\nA2\tA\tIND_A.2\nB1\tB\tIND_B.1\n')
        self.gfa = self.root / 'clip.gfa'
        self.gfa.write_text('S\t1\tA\nW\tREF_A1\t0\tchr1\t0\t1\t>1\nW\tIND_A\t2\tchr1\t0\t1\t>1\nW\tIND_B\t1\tchr1\t0\t1\t>1\n')

    def tearDown(self):
        self.temp.cleanup()

    def test_growth_rejects_successful_but_empty_node_output(self):
        (self.root / 'denominators.json').write_text(json.dumps({'haplotype': 3, 'individual': 2}))
        (self.root / 'haplotype.growth.tsv').write_text('count\tnode\ncoverage\t1\nquorum\t0\n0\t0\n')
        with self.assertRaisesRegex(ValueError, 'count type or cohort range'):
            validate_growth(self.root)

    def test_growth_endpoints(self):
        (self.root / 'denominators.json').write_text(json.dumps({'haplotype': 2, 'individual': 2}))
        for unit in ('haplotype', 'individual'):
            (self.root / (unit + '.hist.tsv')).write_text('count\tbp\n0\t0\n1\t10\n2\t20\n')
            (self.root / (unit + '.growth.tsv')).write_text('count\tbp\tbp\ncoverage\t1\t1\nquorum\t0\t1\n0\t0\t0\n1\t25\t25\n2\t30\t20\n')
        validate_growth(self.root)
        path = self.root / 'individual.growth.tsv'
        path.write_text(path.read_text().replace('2\t30\t20', '2\t31.1\t20'))
        with self.assertRaisesRegex(ValueError, 'endpoint'):
            validate_growth(self.root)

    def test_growth_rejects_zero_only_bp_output(self):
        # Actual failure after adding -c bp: TSV was still parsed as an empty GFA.
        (self.root / 'denominators.json').write_text(json.dumps({'haplotype': 10, 'individual': 5}))
        (self.root / 'haplotype.growth.tsv').write_text(
            'panacus\tgrowth\tgrowth\tgrowth\ncount\tbp\tbp\tbp\n'
            'coverage\t1\t1\t1\nquorum\t0\t1\t0.9\n0\t0\t0\t0\n')
        with self.assertRaisesRegex(ValueError, 'count type or cohort range'):
            validate_growth(self.root)

    def test_reference_and_partner_share_individual(self):
        meta = prepare(self.gfa, self.ledger, self.root)
        self.assertEqual((meta['haplotype'], meta['individual']), (3, 2))
        groups = dict(line.split('\t') for line in (self.root / 'individual.groups.tsv').read_text().splitlines())
        self.assertEqual(groups['REF_A1#0#chr1'], groups['IND_A#2#chr1'])

    def test_unknown_path_fails(self):
        self.gfa.write_text(self.gfa.read_text() + 'W\tunknown\t0\tchr1\t0\t1\t>1\n')
        with self.assertRaisesRegex(ValueError, 'absent from identity'):
            prepare(self.gfa, self.ledger, self.root)

    def test_missing_member_fails(self):
        self.gfa.write_text('W\tREF_A1\t0\tchr1\t0\t1\t>1\n')
        with self.assertRaisesRegex(ValueError, 'missing from graph'):
            prepare(self.gfa, self.ledger, self.root)

    def test_zero_core_does_not_reduce_denominator(self):
        path = self.root / 'hist.tsv'
        path.write_text('panacus\thist\ncount\tbp\n0\t0\n1\t10\n2\t20\n3\t0\n')
        self.assertEqual(read_hist(path, 3)[3], 0)
        with self.assertRaises(ValueError):
            read_hist(path, 2)

    def test_two_individuals_and_overlapping_private(self):
        (self.root / 'denominators.json').write_text(json.dumps({'haplotype': 3, 'individual': 2}))
        (self.root / 'haplotype.hist.tsv').write_text('panacus\thist\ncount\tbp\n0\t0\n1\t10\n2\t20\n3\t0\n')
        (self.root / 'individual.hist.tsv').write_text('panacus\thist\ncount\tbp\n0\t0\n1\t10\n2\t20\n')
        summarize(self.root, 1, -1, 2)
        with (self.root / 'sharing_summary.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(next(r['bp'] for r in rows if r['unit'] == 'haplotype' and r['class'] == 'core'), '0')
        self.assertEqual(next(r['bp'] for r in rows if r['unit'] == 'individual' and r['class'] == 'private'), '10')


if __name__ == '__main__':
    unittest.main()
