import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts/comparisons'))
from sv_alignment_qc import assess, boundary_gap_distance, fasta_index, interval, union_size


class SVAlignmentQCTests(unittest.TestCase):
    def test_reverse_and_overlap(self):
        self.assertEqual(interval(20, 11), (10, 20))
        self.assertEqual(union_size([(10, 20), (15, 25), (25, 30)]), 20)

    def test_gap_across_lines(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'a.fa'
            p.write_text('>x\nAANN\nNnAA\n')
            lengths, gaps, digest = fasta_index(p)
        self.assertEqual(lengths, {'x': 8})
        self.assertEqual(gaps['x'], [(2, 6)])
        self.assertEqual(boundary_gap_distance((3, 7), gaps['x']), 0)
        self.assertEqual(len(digest), 64)

    def test_assessment_preserves_id_and_missing_evidence(self):
        e = dict(ref='x', rstart='101', rend='200', query='x', qstart='300',
                 qend='201', id='INV1', parent='-', type='INV', copy_status='-')
        index = ({'x': 1000}, {'x': []}, '')
        child = dict(e, id='AL1', parent='INV1', type='INVAL')
        result = assess(e, [child, child], index, index, 10, .5)
        self.assertEqual(result['reference_assigned_span_fraction'], 1)
        self.assertEqual(result['qc_status'], 'NO_FLAGS_IN_IMPLEMENTED_CHECKS')
        self.assertEqual(result['read_support'], 'NOT_ASSESSED')
        missing = assess(e, [], index, index, 10, .5)
        self.assertIsNone(missing['reference_assigned_span_fraction'])
        self.assertEqual(missing['id'], 'INV1')
        self.assertIn('unavailable', missing['qc_flags'])

    def test_out_of_bounds_rejected(self):
        e = dict(ref='x', rstart='1', rend='200', query='x', qstart='1',
                 qend='200', id='INV1', parent='-', type='INV', copy_status='-')
        index = ({'x': 100}, {'x': []}, '')
        with self.assertRaises(ValueError):
            assess(e, [], index, index, 10, .5)


if __name__ == '__main__':
    unittest.main()
