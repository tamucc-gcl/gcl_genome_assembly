"""Synthetic regressions for alignment-supported transitions; run remotely."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_scripts"))
from chimera_intervals import regions, intervals, compatible_gaps, verdict_for_interval, arm_context


class IntervalTests(unittest.TestCase):
    def test_out_and_back_is_visible_without_promoting_a_cut(self):
        alns = [(0, 1000, "chr1"), (1010, 1020, "chr2"), (1030, 2030, "chr1")]
        changes = intervals(alns, 1)
        context = arm_context(alns, changes[0], 1)
        self.assertEqual(context["arm_pattern"], "out_and_back")
        self.assertEqual(context["right_arm_aligned_bp"], 10)
        self.assertLess(context["right_arm_fraction"], 0.01)
        self.assertEqual(len(changes), 2)

    def test_arm_support_is_union_not_bounding_span_or_duplicate_sum(self):
        alns = [(0, 100, "chr1"), (0, 100, "chr1"), (200, 300, "chr1"),
                (400, 600, "chr2")]
        context = arm_context(alns, intervals(alns, 1)[0], 1)
        self.assertEqual(context["left_arm_aligned_bp"], 200)
        self.assertEqual(context["left_arm_fraction"], 0.5)
        self.assertEqual(context["arm_pattern"], "chromosome_transition")

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
