"""Remote synthetic tests; these tests never invoke minimap2 or use sequencing data."""
import argparse
import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_scripts"))
from recover_older_joins import prepare, exact_hits, locate, assess, checksum, rows, chromosome


class OlderJoinTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cwd = os.getcwd()
        os.chdir(self.tmp.name)
        Path("original.fa").write_text(">old\nAAAACCCCNNNNGGGGTTTT\n")
        Path("original.agp").write_text(
            "old\t1\t8\t1\tW\tl\t1\t8\t+\n"
            "old\t9\t12\t2\tN\t4\tscaffold\tyes\tproximity_ligation\n"
            "old\t13\t20\t3\tW\tr\t1\t8\t+\n")
        prepare("original.agp", "original.fa", 4, 2, "a")

    def tearDown(self):
        os.chdir(self.cwd)
        self.tmp.cleanup()

    def hit(self, side, start, strand="+", target="new", quality=60):
        return f"J00000000_{side}\t4\t0\t4\t{strand}\t{target}\t22\t{start}\t{start+4}\t4\t4\t{quality}\tcg:Z:4=\n"

    def joint(self):
        return json.loads(Path("a.source.json").read_text())["joins"][0]

    def test_insertion_before_join_moves_coordinate(self):
        Path("flanks.paf").write_text(self.hit("L", 6)+self.hit("R", 14))
        status, loc = locate(self.joint(), exact_hits("flanks.paf"), 30)
        self.assertEqual(status, "mapped_pending_gap_check")
        self.assertEqual(loc, ("new", 10, 14, "+"))

    def test_reverse_orientation(self):
        Path("flanks.paf").write_text(self.hit("L", 14, "-")+self.hit("R", 6, "-"))
        self.assertEqual(locate(self.joint(), exact_hits("flanks.paf"), 30)[1],
                         ("new", 10, 14, "-"))

    def test_duplicate_secondary_rejected_even_low_mapq(self):
        Path("flanks.paf").write_text(self.hit("L", 6)+self.hit("L", 6, target="duplicate", quality=0)+self.hit("R", 14))
        self.assertEqual(locate(self.joint(), exact_hits("flanks.paf"), 30)[0],
                         "unresolved_duplicate_flank")

    def test_split_flanks_unresolved(self):
        Path("flanks.paf").write_text(self.hit("L", 6)+self.hit("R", 14, target="other"))
        self.assertEqual(locate(self.joint(), exact_hits("flanks.paf"), 30)[0],
                         "unresolved_split_or_reoriented_flanks")

    def test_changed_or_missing_flank_unresolved(self):
        Path("flanks.paf").write_text(self.hit("L", 6).replace("cg:Z:4=", "cg:Z:2=1X1="))
        self.assertEqual(locate(self.joint(), exact_hits("flanks.paf"), 30)[0],
                         "unresolved_missing_or_changed_flank")

    def test_source_gap_must_match_fasta(self):
        Path("original.fa").write_text(">old\nAAAACCCCAAAAGGGGTTTT\n")
        with self.assertRaisesRegex(ValueError, "AGP gap"):
            prepare("original.agp", "original.fa", 4, 2, "bad")

    def test_alignment_overlap_not_double_counted(self):
        alns = [(0, 8, "chr1"), (0, 8, "chr1"), (0, 7, "chr2")]
        self.assertEqual(chromosome(alns, 0, 8, 4, 2), ".")

    def assess_fixture(self, gap="NNNN"):
        Path("current.fa").write_text(">new\nTTAAAACCCC"+gap+"GGGGTTTT\n")
        Path("flanks.paf").write_text(self.hit("L", 6)+self.hit("R", 14))
        native = dict(assembly="a", scaffold="new", name="chr1_1+chr2_1", cut_bp="12",
                      left_chrom="chr1", right_chrom="chr2", callable="no", evidence_only="yes",
                      candidate_verdict="REVIEW", chromosome_member="yes", transition_lo="10",
                      transition_hi="14", location_status="unresolved", agp_join_bp=".",
                      agp_join_distance=".", agp_source=".", gap_len=".", reason="no_gap_between_alignment_anchors",
                      join_scope="round2_only", assembly_sha256=checksum("current.fa"),
                      coordinate_stage="pre_finishing")
        with open("native.tsv", "w") as out:
            writer = csv.DictWriter(out, fieldnames=list(native), delimiter="\t")
            writer.writeheader()
            writer.writerow(native)
        Path("map.tsv").write_text("old_name\tnew_name\nr1\tchr1_1\nr2\tchr2_1\n")
        Path("ref.paf").write_text(
            "new\t22\t0\t10\t+\tr1\t50\t0\t10\t10\t10\t60\n"
            "new\t22\t14\t22\t+\tr2\t50\t0\t8\t8\t8\t60\n")
        Path("candidates.tsv").write_text("assembly\tscaffold\tname\tchromosome_member\tverdict\na\tnew\tchr1_1+chr2_1\tyes\tREVIEW\n")
        assess(argparse.Namespace(prefix="a", flank_paf="flanks.paf", ref_map="map.tsv",
            ref_paf="ref.paf", candidates="candidates.tsv", assembly="a", min_mapq=30,
            current="current.fa", native="native.tsv", window=8, min_bp=4, margin=2,
            min_block=1, min_span=0, transition_max_bridge=10))

    def test_recovered_join_is_review_only_and_stamped(self):
        self.assess_fixture()
        call = rows("a.review_joins.tsv")[1][0]
        self.assertEqual(call["cut_bp"], "12")
        self.assertEqual(call["candidate_verdict"], "REVIEW")
        self.assertEqual(call["assembly_sha256"], checksum("current.fa"))
        self.assertEqual(rows("a.older_join_audit.tsv")[1][0]["auto_cut"], "no")


    def test_interval_lists_recovered_joins_without_authorizing_cut(self):
        self.assess_fixture()
        transitions = rows("a.transition_intervals.tsv")[1]
        self.assertEqual(transitions[0]["recovered_join_ids"], "J00000000")
        self.assertEqual(transitions[0]["auto_cut"], "no")

    def test_replaced_gap_not_promoted_to_safe_join(self):
        self.assess_fixture("ACGT")
        calls = rows("a.review_joins.tsv")[1]
        self.assertTrue(calls)
        self.assertTrue(all(r["callable"] == "no" and r["evidence_only"] == "yes" for r in calls))
        self.assertEqual(rows("a.older_join_audit.tsv")[1][0]["status"],
                         "unresolved_gap_changed_or_replaced")


if __name__ == "__main__":
    unittest.main()
