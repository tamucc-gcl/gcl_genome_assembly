import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts/comparisons'))
from inspect_gref_catalog import choose_windows


class GrefWindowTests(unittest.TestCase):
    def test_short_and_overlapping_windows_merge(self):
        rows = [('chr1_1', 15000, 5), ('chr1_1_1_alt', 100, 2)]
        self.assertEqual(choose_windows(rows), [('chr1_1', 0, 15000, 'reference_named'),
                                              ('chr1_1_1_alt', 0, 100, 'alt_named')])

    def test_long_coordinates_are_bounded_and_deterministic(self):
        rows = [(f'chr{i}_1', 1000000, 10) for i in range(1, 101)]
        result = choose_windows(rows)
        self.assertEqual(result, choose_windows(list(reversed(rows))))
        self.assertEqual(len(result), 36)
        self.assertTrue(all(end - start == 10000 for _, start, end, _ in result))


if __name__ == '__main__':
    unittest.main()
