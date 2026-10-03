"""Synthetic support-scan regressions; execute with the remote test suite."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('support_scan', Path(__file__).resolve().parents[1] /
                                             'scripts/comparisons/scan_chimera_support.py')
scan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan)


class SupportScanTests(unittest.TestCase):
    def record(self, cigar, mq=60, nm=1, flag=0):
        return ['read', str(flag), 'scaffold', '1', str(mq), cigar, '*', '0', '0', '*', '*', 'NM:i:'+str(nm)]

    def test_small_indel_does_not_destroy_chain(self):
        start, end, blocks, rate = scan.alignment(self.record('1500M1I1500M'))
        self.assertEqual((start, end, blocks), (0, 3000, [(0, 3000)]))
        self.assertAlmostEqual(rate, 1/3001)

    def test_large_deletion_is_not_continuity(self):
        self.assertEqual(scan.alignment(self.record('1500M100D1500M'))[2],
                         [(0, 1500), (1600, 3100)])

    def test_large_insertion_splits_chain(self):
        self.assertEqual(scan.alignment(self.record('1500M100I1500M'))[2],
                         [(0, 1500), (1500, 3000)])

    def test_low_mapq_support_retained_but_not_promoted(self):
        with tempfile.TemporaryDirectory() as tmp:
            sam = Path(tmp)/'a.sam'
            sam.write_text('\t'.join(self.record('4000M', mq=10))+'\n')
            names, rows = scan.scan(sam, dict(scaffold='scaffold', start=0, end=4000,
                                             lo=1000, hi=3000), 'A'*4000)
            self.assertEqual(names, {'read'})
            self.assertEqual(rows[2]['primary_span_any_mapq'], 1)
            self.assertEqual(rows[2]['primary_span_mapq20'], 0)
            self.assertEqual(rows[2]['chain_mapq20_nm1pct'], 0)

    def test_high_mapq_high_error_is_not_clean_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            sam = Path(tmp)/'a.sam'
            sam.write_text('\t'.join(self.record('4000M', nm=800))+'\n')
            _, rows = scan.scan(sam, dict(scaffold='scaffold', start=0, end=4000,
                                         lo=1000, hi=3000), 'A'*4000)
            self.assertEqual(rows[2]['primary_span_mapq20'], 1)
            self.assertEqual(rows[2]['chain_mapq20_nm2pct'], 0)


if __name__ == '__main__':
    unittest.main()
