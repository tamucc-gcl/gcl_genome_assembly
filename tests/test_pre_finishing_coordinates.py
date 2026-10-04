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

    def test_empty_reference_and_populated_tables_have_identical_stamped_headers(self):
        self.guard()
        populated_header = self.table.read_text().splitlines()[0]
        self.run_script("chimera_schema.py", "--empty", self.table)
        self.run_script("chimera_coordinate_guard.py", "--fasta", self.fa,
            "--agp", self.agp, "--table", self.table, "--audit", self.d / "audit.tsv",
            "--round", "round2", "--alignment-status", "unavailable")
        self.assertEqual(self.table.read_text().splitlines()[0], populated_header)
        self.assertIn("alignment_transition_assessment\tunavailable",
                      (self.d / "audit.tsv").read_text())
        self.guard()
        fields = self.table.read_text().splitlines()[0].split("\t")
        self.assertEqual(len(fields), len(set(fields)))

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

    def cut_fixture(self, member="yes", diagnostic="no", mode="file", ok=True):
        self.table.write_text(
            "assembly\tscaffold\tcut_bp\tleft_chrom\tright_chrom\tcallable\t"
            "candidate_verdict\tchromosome_member\tevidence_only\tlocation_status\t"
            "transition_lo\ttransition_hi\tagp_join_bp\n"
            f"asm\ts1\t13\tchr1\tchr2\tyes\tBREAK_CANDIDATE\t{member}\t{diagnostic}\t"
            "unique_supported_gap\t10\t15\t13\n")
        self.guard()
        nm = self.d / "names.tsv"
        nm.write_text("old_name\tnew_name\tlength\tclass\tchromosome_member\n"
                      f"s1\tchr1_1+chr2_1\t25\tcomposite\t{member}\n"
                      "s2\tchr3_1\t10\tchromosome\tyes\n")
        return self.run_script("break_chimeras.py", "--fasta", self.fa,
            "--name-map", nm, "--candidates", self.table, "--assembly", "asm",
            "--out-fasta", self.d / "out.fa", "--out-name-map", self.d / "out.tsv",
            "--mode", mode, "--min-piece-bp", "1", ok=ok)

    def test_supplied_file_cannot_bypass_chromosome_scope(self):
        result = self.cut_fixture(member="no", ok=False)
        self.assertIn("Unsafe supplied breakpoint", result.stderr)

    def test_supplied_file_cannot_bypass_unresolved_inference(self):
        self.cut_fixture(member="unresolved", ok=False)

    def test_diagnostic_position_cannot_be_supplied_as_cut(self):
        self.cut_fixture(diagnostic="yes", ok=False)

    def test_auto_rejects_even_a_scoped_vote_candidate(self):
        result = self.cut_fixture(mode="auto", ok=False)
        self.assertIn("private/heterozygous", result.stderr)
        self.assertFalse((self.d / "out.fa").exists())

    def test_cross_scaffold_partner_is_preserved_with_lifted_orientation(self):
        self.pairs.write_text("r3\tA\t12\tA\t22\t+\t-\tUU\n")
        self.lift()
        import gzip
        with gzip.open(self.d / "asm.alternative_partners.pairs.gz", "rt") as handle:
            text = handle.read()
        self.assertIn("r3\ts1\t24\ts2\t2\t-\t-\tUU", text)
        self.assertEqual((self.d / "asm.s1.pairs").read_text(), "")
        self.assertIn("alternative_partner_pairs_written\t1",
                      (self.d / "asm.hic_pairs_audit.tsv").read_text())

    def test_valid_scoped_gap_can_be_cut(self):
        self.cut_fixture()
        self.assertIn(">s1_sub_0_13", (self.d / "out.fa").read_text())


if __name__ == "__main__":
    unittest.main()
