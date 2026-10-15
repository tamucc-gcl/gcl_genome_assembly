"""Report stage ordering and metric fidelity, independent of study assemblies."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_scripts"))
from generate_assembly_report import parse_metrics, get_available_stages, build_report_data


class AssemblyQCReportTests(unittest.TestCase):
    def test_pipeline_stage_order_and_disabled_stages(self):
        expected = ["ctg.base", "ctg.org", "ctg.purged", "ctg.cor", "ctg.deco",
                    "scaf.base", "scaf.cor", "scaf.deco", "scaf2", "chim.brk",
                    "gap_fill", "teloclip", "final"]
        data = {"sample": {stage: {} for stage in reversed(expected)}}
        self.assertEqual(get_available_stages(data), expected)
        self.assertEqual(get_available_stages({"sample": {"final": {}}}), ["final"])

    def test_haplotypes_and_combined_completeness_are_distinct(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "metrics.csv"
            with path.open("w", newline="", encoding="utf-8") as dest:
                writer = csv.writer(dest)
                writer.writerow(["sample_id", "stage", "analysis", "metric", "hap1", "hap2", "both"])
                writer.writerows([
                    ["arbitrary", "final", "quast", "N50", 8000000, 9000000, "NA"],
                    ["arbitrary", "final", "busco", "complete", 95, 90, "NA"],
                    ["arbitrary", "final", "busco", "total_busco", 100, 100, "NA"],
                    ["arbitrary", "final", "merqury", "kmer_completeness", 76, 78, 98],
                ])
            data = parse_metrics(path)
            report = build_report_data(data, get_available_stages(data))
            assembly = report["stage_data"]["final"]["assemblies"]["arbitrary"]
            self.assertEqual(assembly["hap1"]["n50"], 8)
            self.assertEqual(assembly["hap2"]["busco_pct"], 90)
            self.assertEqual(assembly["both"]["kcomp"], 98)
            self.assertNotIn("tc_ratio", assembly["hap1"])


if __name__ == "__main__":
    unittest.main()
