"""Inference-to-chimera-scope regressions; synthetic lengths only."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_scripts"))
from harmonize_names import select_chromosome_set, chromosome_scope_status


class ChromosomeScopeTests(unittest.TestCase):
    def infer(self, lengths, method="dropoff"):
        selected, meta = select_chromosome_set(
            lengths, min_scaffold_bp=1, method=method, dropoff_ratio=2.0,
            dropoff_min_frac=0.5, min_chrom_frac=0.0)
        return {name for name, _ in selected}, meta

    def test_failed_dropoff_keeps_naming_set_but_not_cut_scope(self):
        selected, meta = self.infer([("a", 100), ("b", 95), ("c", 90), ("d", 85)])
        self.assertEqual(meta["method"], "threshold_fallback")
        self.assertEqual(len(selected), 4)
        for scaffold in selected:
            member, reason = chromosome_scope_status(scaffold, selected, meta)
            self.assertEqual(member, "unresolved")
            self.assertIn("threshold_fallback", reason)

    def test_successful_dropoff_keeps_only_inferred_members(self):
        selected, meta = self.infer([("a", 100), ("b", 90), ("tail", 2)])
        self.assertEqual(meta["method"], "dropoff")
        self.assertEqual(chromosome_scope_status("a", selected, meta)[0], "yes")
        self.assertEqual(chromosome_scope_status("tail", selected, meta)[0], "no")

    def test_nonvoter_with_inferred_chromosomes_remains_eligible(self):
        selected, meta = self.infer([("x", 100), ("y", 90), ("z", 2)])
        meta["role"] = "passenger"
        self.assertEqual(chromosome_scope_status("x", selected, meta)[0], "yes")

    def test_threshold_only_selection_is_not_inference(self):
        selected, meta = self.infer([("a", 100), ("b", 90)], method="threshold")
        self.assertEqual(chromosome_scope_status("a", selected, meta)[0], "unresolved")

    def test_missing_quality_metadata_fails_closed(self):
        self.assertEqual(chromosome_scope_status("a", {"a"}, {})[0], "unresolved")

    def test_scope_does_not_depend_on_absolute_length_or_sample_name(self):
        for scale in (1, 1000000):
            selected, meta = self.infer([("renamed_X", 100*scale), ("renamed_Y", 90*scale), ("tail", 2*scale)])
            self.assertEqual(chromosome_scope_status("renamed_X", selected, meta)[0], "yes")


if __name__ == "__main__":
    unittest.main()
