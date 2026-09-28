"""Dependency-free control-flow tests; numerical/plot libraries are mocked."""
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

SCRIPT = Path(__file__).resolve().parents[1] / "py_scripts" / "chimera_evidence.py"


class EvidenceReportingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        spec = importlib.util.spec_from_file_location("evidence_under_test", SCRIPT)
        self.module = importlib.util.module_from_spec(spec)
        # The cluster's base Python need not have numpy/matplotlib/cooler installed.
        with patch.dict(sys.modules, {"numpy": Mock()}):
            spec.loader.exec_module(self.module)
        self.row = dict(assembly="a", scaffold="s", cut_bp="50",
                        candidate_verdict="REVIEW", name="s", vote="not_recomputed",
                        evidence_only="no")
        self.cool = self.root / "x.cool"
        self.cool.touch()

    def run_main(self, diagnostic=False, make_figure=True, requested=50):
        m = self.module
        row = dict(self.row)
        if diagnostic:
            row.update(evidence_only="yes", transition_lo="20", transition_hi="80")
        hic = dict(min_bp=40, ratio=0.5, min_value=1, median=2,
                   n_low_contiguous=2, low_bp=[40, 50])
        def draw(path, *args, **kwargs):
            if make_figure:
                Path(path).write_bytes(b"mock figure")
        args = ["chimera_evidence.py", "--fasta", "unused.fa", "--candidates", "unused.tsv",
                "--assembly", "a", "--scaffold", "s", "--cut-bp", str(requested),
                "--outdir", str(self.root), "--label", "a.s_50", "--cool", str(self.cool)]
        with patch.object(sys, "argv", args), \
             patch.object(m, "read_rows", return_value=[row]), \
             patch.object(m, "read_one_fasta", return_value=("s", "A"*100)), \
             patch.object(m, "hic_profile", return_value=([1, 2], 10, "")), \
             patch.object(m, "summarise_profile", return_value=hic), \
             patch.object(m, "ratio_at", return_value=dict(ratio=0.8, value=1.6, note="")), \
             patch.object(m, "telomere_windows", return_value=[]), \
             patch.object(m, "summarise_telomere", return_value=None), \
             patch.object(m, "gaps", return_value=[]), \
             patch.object(m, "snap", return_value=None) as snap_mock, \
             patch.object(m, "figure", side_effect=draw):
            m.main()
            if diagnostic:
                snap_mock.assert_not_called()
        text = (self.root / "a.s_50.chimera_evidence.tsv").read_text()
        return text

    def test_metrics_loop_preserves_figure_prefix_and_verdict(self):
        text = self.run_main()
        self.assertTrue((self.root / "a.s_50.chimera_evidence.png").is_file())
        self.assertFalse((self.root / "n_low_contiguous.chimera_evidence.png").exists())
        self.assertIn("verdict_from_candidates\tREVIEW", text)
        self.assertIn("figure_status\tgenerated", text)

    def test_diagnostic_position_is_not_snapped_or_reported_as_cut(self):
        text = self.run_main(diagnostic=True)
        self.assertIn("cut_bp\tNA", text)
        self.assertIn("evidence_position_bp\t50", text)
        self.assertIn("transition_lo\t20\ntransition_hi\t80", text)
        self.assertNotIn("no_gap_within", text)

    def test_absent_requested_position_fails_instead_of_using_first_row(self):
        with self.assertRaises(SystemExit):
            self.run_main(requested=51)

    def test_missing_expected_figure_fails(self):
        with self.assertRaisesRegex(RuntimeError, "Expected evidence figure"):
            self.run_main(make_figure=False)


if __name__ == "__main__":
    unittest.main()
