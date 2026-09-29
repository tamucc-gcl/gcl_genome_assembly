"""Exercise the caller with tiny, sample-neutral files; no aligner is invoked."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_scripts"))
import chimera_joins


class CandidateScopeTests(unittest.TestCase):
    def call(self, member="yes", gaps=(55,), verdict="BREAK_CANDIDATE"):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def write(name, text):
                p = root / name
                p.write_text(text)
                return str(p)
            candidates = write("candidates.tsv",
                "assembly\tscaffold\tname\tspan_bp\tmembers\tvote\tverdict\tchromosome_member\n"
                f"a\ts\tchr1_1+chr2_1\t100\tref1+ref2\t1f/3s\t{verdict}\t{member}\n")
            names = write("names.tsv", "old_name\tnew_name\nr1\tchr1_1\nr2\tchr2_1\n")
            paf = write("ref.paf", "s\t100\t0\t50\t+\tr1\t100\t0\t50\t50\t50\t60\n"
                        "s\t100\t60\t100\t-\tr2\t100\t0\t40\t40\t40\t60\n")
            joins = write("joins.tsv", "assembly\tfinal_object\tfinal_cut\tgap_len\tsource\n" +
                          "".join(f"a\ts\t{pos}\t10\tcurrent\n" for pos in gaps))
            output = root / "out.tsv"
            args = ["chimera_joins.py", "--assembly", "a", "--scaffold-map", candidates,
                    "--ref-name-map", names, "--paf", paf, "--joins", joins,
                    "--round1", "unused.agp", "--min-block", "1", "--out", str(output)]
            with patch.object(sys, "argv", args):
                chimera_joins.main()
            with output.open() as handle:
                return list(csv.DictReader(handle, delimiter="\t"))

    def test_small_inferred_chromosome_is_eligible(self):
        row = self.call()[0]
        self.assertEqual(row["callable"], "yes")
        self.assertEqual(row["cut_bp"], "55")
        self.assertEqual(row["candidate_verdict"], "BREAK_CANDIDATE")

    def test_excluded_scaffold_produces_no_profile(self):
        self.assertEqual(self.call(member="no"), [])

    def test_non_candidate_produces_no_profile(self):
        self.assertEqual(self.call(verdict="NOT_A_CANDIDATE"), [])

    def test_upstream_gap_is_not_snapped_into_transition(self):
        row = self.call(gaps=(35,))[0]
        self.assertEqual(row["callable"], "no")
        self.assertEqual(row["evidence_only"], "yes")
        self.assertEqual(row["cut_bp"], "55")

    def test_multiple_gaps_require_review(self):
        row = self.call(gaps=(52, 58))[0]
        self.assertEqual(row["candidate_verdict"], "REVIEW")
        self.assertEqual(row["callable"], "no")


if __name__ == "__main__":
    unittest.main()
