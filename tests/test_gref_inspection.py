import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts/comparisons'))
from inspect_gref_catalog import choose_windows, ledger_slots, slot_accounting


class GrefWindowTests(unittest.TestCase):
    def test_reference_partner_slots_are_not_extra_haplotypes(self):
        slots = ledger_slots([{'graph_name': 'REF_reference'}, {'graph_name': 'IND_partner.2'}])
        self.assertEqual(slots, {'REF_reference': {0}, 'IND_partner': {1}})
        result = slot_accounting('.|1', slots['IND_partner'])
        self.assertEqual(result['biological_called'], 1)
        self.assertEqual(result['placeholder_missing'], 1)
        self.assertEqual(result['biological_missing'], 0)

    def test_called_placeholder_is_visible(self):
        self.assertEqual(slot_accounting('0|1', {1})['placeholder_called'], 1)

    def test_unphased_partner_and_short_gt_are_unresolved(self):
        for gt in ('./1', '1', '.'):
            self.assertEqual(slot_accounting(gt, {1})['unresolved_records'], 1)

    def test_ordinary_diploid_missingness(self):
        counts = slot_accounting('0/.', {0, 1})
        self.assertEqual(counts['biological_called'], 1)
        self.assertEqual(counts['biological_missing'], 1)
        self.assertEqual(counts['placeholder_missing'], 0)

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
