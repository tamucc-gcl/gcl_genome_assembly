import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/comparisons'))
from compare_chimera_history import summarize, union_length


class HistorySummaryTests(unittest.TestCase):
    def test_union_does_not_double_count_secondary_overlaps(self):
        self.assertEqual(union_length([(10, 20), (0, 15), (30, 40), (3, 8)]), 30)

    def test_partners_and_orientations_remain_separate(self):
        base = dict(query='new', target='old1', strand='+', query_start=0,
                    query_end=100, target_start=200, target_end=300,
                    query_length=1000, target_length=2000)
        rows = summarize([base, dict(base), dict(base, target='old2'), dict(base, strand='-')])
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(r['query_covered_bp'] == 100 for r in rows))
        self.assertEqual(sum(r['blocks'] for r in rows), 4)


if __name__ == '__main__':
    unittest.main()
