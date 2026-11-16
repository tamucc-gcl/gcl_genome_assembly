import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_scripts"))
from contamination_evidence import (blastn_result, followup_plan, locus_comparison, ranked_loci,
                                    search_hits, whole_contig_evidence)


def hit(start, end, relation="non-host", group="Arthropoda"):
    return dict(start0=start, end0=end, host_relation=relation,
                lineage_at_host_rank=group, bitscore=100, subject="accession", identity_pct=92)


class ContaminationEvidenceTests(unittest.TestCase):
    def test_short_or_competing_hits_stay_visible_at_lower_priority(self):
        loci = ranked_loci([hit(100, 200), hit(600, 1900), hit(600, 1900, "host-compatible")], [])
        self.assertEqual([locus["priority"] for locus in loci],
                         ["short non-host signal", "competing hits only"])
        self.assertEqual([locus["distinct_bp"] for locus in loci], [100, 0])

    def test_concordance_requires_matching_taxonomy_and_overlapping_hits(self):
        protein = [hit(200, 700)]
        same = ranked_loci(protein, [hit(300, 800)])[0]
        other = ranked_loci(protein, [hit(300, 800, group="Mollusca")])[0]
        self.assertEqual(same["priority"], "concordant non-host")
        self.assertEqual(other["priority"], "short non-host signal")
        self.assertEqual(same["distinct_bp"], 600)

    def test_followup_budget_keeps_cropped_concordant_locus_first(self):
        windows = {
            "extended": dict(sequence="a", start0="0", end0="100000"),
            "concordant": dict(sequence="b", start0="100000", end0="200000"),
            "short": dict(sequence="c", start0="0", end0="100000"),
        }
        protein = {"extended": [hit(5000, 6500)], "concordant": [hit(100100, 100700)],
                   "short": [hit(8000, 8200)]}
        mega = {"concordant": [hit(100150, 100750)]}
        plan = followup_plan(windows, protein, mega, 1250)
        self.assertEqual(list(plan), ["concordant"])
        self.assertEqual((plan["concordant"]["start0"], plan["concordant"]["end0"]),
                         (100000, 101250))

    def test_same_database_agreement_does_not_promote_locus(self):
        locus, = ranked_loci([], [hit(200, 700)], [hit(300, 800)])
        self.assertEqual(locus["priority"], "short non-host signal")
        self.assertIn("same nt database", locus_comparison(locus, {})["search_agreement"])

    def test_comparison_reports_the_host_hit_on_the_same_bases(self):
        foreign = hit(100, 1500)
        local = hit(300, 1300, "host-compatible", "Chordata")
        distant = dict(hit(1600, 2500, "host-compatible", "Chordata"), bitscore=400)
        locus, = ranked_loci([], [foreign])
        comparison = locus_comparison(locus, {"megablast": [local, distant]})
        self.assertIn("1,400 bp in interval", comparison["non_host_match"])
        self.assertIn("1,000 bp overlapping non-host hit", comparison["host_match"])
        self.assertNotIn("400 bits", comparison["host_match"])

    def test_whole_contig_screen_separates_interval_hits_from_whole_contig_hits(self):
        spans = {"Arthropoda": [(0, 600), (500, 700)]}
        self.assertEqual(whole_contig_evidence(spans, 700, 50, 1000)["whole_contig_candidate"], "YES")
        self.assertEqual(whole_contig_evidence(spans, 700, 50, 10000)["whole_contig_candidate"], "NO")
        self.assertEqual(whole_contig_evidence(spans, 700, 150, 1000)["whole_contig_candidate"], "NO")
        self.assertEqual(whole_contig_evidence(spans, 900, 50, 1000)["whole_contig_candidate"], "NO")

    def test_followup_can_report_a_host_only_hit_on_part_of_a_nominated_locus(self):
        locus, = ranked_loci([], [hit(100, 1500)])
        target = dict(start0=50, end0=1000, hit_start0=100, hit_end0=1000)
        host = hit(200, 700, "host-compatible", "Chordata")
        status, nominated = blastn_result(locus, target, [host], searched=True)
        self.assertTrue(nominated)
        self.assertEqual(status, "BLASTN: overlapping host hit found (900 bp of interval searched)")

    def test_followup_alignment_uses_cropped_origin_and_preserves_strand(self):
        target = {"query": dict(sequence="contig", start0=1500, end0=2800)}
        args = SimpleNamespace(taxdump=None, host_rank="phylum", taxid=1)
        with tempfile.TemporaryDirectory() as work:
            raw = Path(work) / "blastn.tsv"
            raw.write_text("query\tsubject\t42\t600\t101\t100\t1e-30\t95\ttitle\n")
            with patch("contamination_evidence.taxonomy", return_value=(1, {42: 2}, {2: "Arthropoda"})):
                result, = search_hits([raw], target, args)
        self.assertEqual((result["sequence"], result["start0"], result["end0"], result["strand"]),
                         ("contig", 1600, 2100, "-"))


if __name__ == "__main__":
    unittest.main()
