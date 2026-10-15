"""General scaffold eligibility and sequence preservation checks."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_scripts"))
from teloclip_scope import prepare, verify


class TeloclipScopeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root / "assembly.fa"
        self.fasta.write_text(">large\nAAAA\nCCCC\n>small\nGGGG\n>uncertain\nTTTT\n")
        self.mapping = self.root / "names.tsv"
        self.mapping.write_text("old_name\tchromosome_member\nlarge\tyes\nsmall\tno\nuncertain\tunresolved\n")
        self.prefix = self.root / "audit"

    def test_only_explicit_members_are_eligible_and_nonmembers_preserved(self):
        self.assertEqual(prepare(self.fasta, self.mapping, self.prefix), 1)
        self.assertEqual((self.root / "audit.eligible.bed").read_text(), "large\t0\t8\n")
        extended = self.root / "extended.fa"
        extended.write_text(">large\nTTAAAACCCCTT\n>small\nGGGG\n>uncertain\nTTTT\n")
        verify(extended, self.prefix)
        with (self.root / "audit.teloclip_scope.tsv").open() as source:
            rows = list(csv.DictReader(source, delimiter="\t"))
        self.assertEqual([r["added_bp"] for r in rows], ["4", "0", "0"])

    def test_ineligible_change_or_record_loss_is_rejected(self):
        prepare(self.fasta, self.mapping, self.prefix)
        for text in [">large\nAAAACCCC\n>small\nGGGGTT\n>uncertain\nTTTT\n",
                     ">large\nAAAACCCC\n>uncertain\nTTTT\n"]:
            extended = self.root / "bad.fa"
            extended.write_text(text)
            with self.assertRaises(ValueError):
                verify(extended, self.prefix)

    def test_no_chromosome_or_no_classification_is_pass_through(self):
        self.mapping.write_text(self.mapping.read_text().replace("yes", "no"))
        self.assertEqual(prepare(self.fasta, self.mapping, self.prefix), 0)
        verify(self.fasta, self.prefix)
        sentinel = self.root / "NO_HARMONIZE"
        sentinel.touch()
        self.assertEqual(prepare(self.fasta, sentinel, self.prefix), 0)
        verify(self.fasta, self.prefix)

    def test_stale_or_malformed_map_fails_instead_of_guessing_from_names(self):
        for text in ["old_name\tnew_name\nlarge\tchr1\n",
                     "old_name\tchromosome_member\nlarge\tyes\n",
                     self.mapping.read_text().replace("unresolved", "invalid")]:
            self.mapping.write_text(text)
            with self.assertRaises(ValueError):
                prepare(self.fasta, self.mapping, self.prefix)


if __name__ == "__main__":
    unittest.main()
