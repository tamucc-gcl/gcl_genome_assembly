"""Synthetic regressions for alignment-supported transitions; run remotely."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_scripts"))
from chimera_intervals import regions, intervals, compatible_gaps, verdict_for_interval


class IntervalTests(unittest.TestCase):
    def test_duplicate_alignments_do_not_add_support(self):
        self.assertEqual(regions([(0, 8, "chr1")]*5, 10), [])

    def test_conflicting_overlap_is_unassigned(self):
        self.assertEqual(regions([(0, 20, "chr1"), (5, 15, "chr2")], 1),
                         [(0, 5, "chr1"), (15, 20, "chr1")])

    def test_gap_inside_same_chromosome_block_is_not_a_break(self):
        # Scaled regression: the old component boundary precedes the true transition.
        alns = [(29987, 30623, "chr7"), (30623, 31397, "chr7"),
                (32161, 33665, "chr7"), (33826, 35788, "chr12")]
        transition = intervals(alns, 100)[0]
        self.assertEqual((transition["lo"], transition["hi"]), (33665, 33826))
        self.assertEqual(compatible_gaps(transition, [(30202, {}), (33004, {})]), [])

    def test_multiple_gaps_preserve_ambiguity(self):
        transition = intervals([(0, 10, "chr1"), (20, 30, "chr2")], 5)[0]
        self.assertEqual(len(compatible_gaps(transition, [(12, {}), (18, {})])), 2)

    def test_multijunction_cannot_inherit_scaffold_cut_verdict(self):
        candidate = dict(members="ref1+ref2", name="chr1_1+chr2_1", verdict="BREAK_CANDIDATE")
        interval = dict(left="chr1", right="chr2")
        self.assertEqual(verdict_for_interval(candidate, interval, 2), "REVIEW")
        self.assertEqual(verdict_for_interval(candidate, interval, 1), "BREAK_CANDIDATE")
        self.assertEqual(verdict_for_interval(candidate, dict(left="chr1", right="chr3"), 1), "REVIEW")

    def test_small_chromosome_has_no_absolute_size_gate(self):
        self.assertEqual(len(intervals([(0, 50, "chr1"), (60, 100, "chr2")], 1)), 1)


if __name__ == "__main__":
    unittest.main()
