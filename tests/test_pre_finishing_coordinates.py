"""Small synthetic coordinate checks. Run on the user's machine; no biological inputs."""
import csv
import hashlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class Coordinates(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)
        self.agp = self.d / "last.agp"
        self.agp.write_text(
            "s1\t1\t10\t1\tW\tA\t1\t10\t+\n"
            "s1\t11\t15\t2\tN\t5\tscaffold\tyes\tproximity_ligation\n"
            "s1\t16\t25\t3\tW\tA\t11\t20\t-\n"
            "s2\t1\t10\t1\tW\tA\t21\t30\t+\n")
        self.fa = self.d / "assessment.fa"
        self.fa.write_text(">s1\n" + "A"*10 + "N"*5 + "C"*10 + "\n>s2\n" + "G"*10 + "\n")
        self.pairs = self.d / "input.pairs"
        self.pairs.write_text("r1\tA\t2\tA\t12\t+\t+\tUU\nr2\tA\t22\tA\t29\t+\t-\tUU\n")
        self.table = self.d / "called.tsv"
        self.table.write_text("assembly\tscaffold\tcut_bp\tleft_chrom\tright_chrom\tcallable\n"
                              "asm\ts1\t13\tchr1\tchr2\tyes\n")

    def tearDown(self):
        self.tmp.cleanup()

    def run_script(self, name, *args, ok=True):
        result = subprocess.run([sys.executable, str(ROOT / "py_scripts" / name),
                                 *map(str, args)], capture_output=True, text=True, cwd=self.d)
        if ok:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def lift(self):
        return self.run_script("chimera_hic_pairs.py", "--pairs", self.pairs,
            "--round1", self.agp, "--scaffolds", "s1,s2", "--outdir", self.d, "--label", "asm")

    def guard(self, ok=True):
        return self.run_script("chimera_coordinate_guard.py", "--fasta", self.fa,
            "--agp", self.agp, "--table", self.table, "--audit", self.d / "audit.tsv",
            "--round", "round2", ok=ok)

    def test_split_component_and_reverse_strand(self):
        self.lift()
        self.assertEqual((self.d / "asm.s1.pairs").read_text(),
                         "r1\ts1\t2\ts1\t24\t+\t-\tUU\n")
        self.assertEqual((self.d / "asm.s2.pairs").read_text(),
                         "r2\ts2\t2\ts2\t9\t+\t-\tUU\n")

    def test_duplicate_placement_is_not_assigned_arbitrarily(self):
        with self.agp.open("a") as handle:
            handle.write("s3\t1\t10\t1\tW\tA\t1\t10\t+\n")
        self.lift()
        self.assertEqual((self.d / "asm.s1.pairs").read_text(), "")
        self.assertIn("ambiguous_ends\t1", (self.d / "asm.hic_pairs_audit.tsv").read_text())

    def test_guard_stamps_exact_fasta(self):
        self.guard()
        with self.table.open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["assembly_sha256"], hashlib.sha256(self.fa.read_bytes()).hexdigest())
        self.assertEqual(row["coordinate_stage"], "pre_finishing")
        self.assertEqual(row["join_scope"], "round2_only")

    def test_guard_rejects_gap_replaced_without_length_change(self):
        self.fa.write_text(self.fa.read_text().replace("NNNNN", "AAAAA"))
        self.guard(ok=False)

    def test_reviewed_break_table_rejects_changed_fasta(self):
        self.guard()
        self.fa.write_text(self.fa.read_text().replace("AAAAAAAAAA", "TAAAAAAAAA"))
        nm = self.d / "names.tsv"
        nm.write_text("old_name\tnew_name\torient\torder\tlength\tclass\tref_span\tflags\n"
                      "s1\tchr1_1\tfwd\t1\t25\tchromosome\t.\t.\n"
                      "s2\tchr2_1\tfwd\t2\t10\tchromosome\t.\t.\n")
        result = self.run_script("break_chimeras.py", "--fasta", self.fa,
            "--name-map", nm, "--candidates", self.table, "--assembly", "asm",
            "--out-fasta", self.d / "out.fa", "--out-name-map", self.d / "out.tsv",
            "--mode", "file", "--min-piece-bp", "1", ok=False)
        self.assertIn("does not match", result.stderr)

    def test_empty_support_is_reported(self):
        self.pairs.write_text("")
        self.lift()
        self.assertIn("no_supporting_pairs", (self.d / "asm.hic_pairs_audit.tsv").read_text())


if __name__ == "__main__":
    unittest.main()