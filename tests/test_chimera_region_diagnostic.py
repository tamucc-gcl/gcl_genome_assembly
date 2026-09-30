"""Synthetic tests of exploratory bin-pair accounting; no sequencing files."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/comparisons'))
from explore_chimera_region import local_contacts


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
