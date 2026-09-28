"""Synthetic eligibility checks for the remote host; no sequencing data required."""
import copy
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_scripts"))
from assembly_eligibility import evaluate

SETTINGS = dict(min_individuals=2, run_pangenome=False, qc_mode="none",
                min_scaffold_bp=1000000, chromosome_method="dropoff",
                dropoff_ratio=2.0, dropoff_min_frac=0.5, min_chrom_frac=0.1,
                min_cut_ratio=0.0, min_genome_frac=0.8, min_n50_ratio=0.2,
                require_dropoff=True, min_voters=3)


class EligibilityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.previous = os.getcwd()
        os.chdir(self.tmp.name)
        self.data = dict(settings=copy.deepcopy(SETTINGS), expected=[], observed=[],
                         harmonization_references=[])

    def tearDown(self):
        os.chdir(self.previous)
        self.tmp.cleanup()

    def add(self, sample, taxid="1", ploidy=2, assembler="hifiasm",
            phased=True, fragmented=False, missing=()):
        n_hap = 2 if assembler == "hifiasm" and ploidy == 2 else 1
        self.data["expected"].append(dict(sample=sample, taxid=taxid, species="Species",
            assembler=assembler, ploidy=ploidy, n_hap=n_hap, hic_phasing=phased))
        for hap in (["hap1", "hap2"] if n_hap == 2 else ["primary"]):
            if hap in missing:
                continue
            aid = sample + "_" + hap
            lengths = [100000] * 30 if fragmented else [10000000, 9000000, 100000]
            fai = aid + ".fasta.fai"
            Path(fai).write_text("".join(f"seq{i}\t{n}\t0\t60\t61\n"
                                         for i, n in enumerate(lengths)))
            self.data["observed"].append(dict(meta=dict(id=aid, sample=sample, taxid=taxid),
                fasta=str(Path(aid + ".fasta").absolute()), fai=str(Path(fai).absolute()),
                fai_name=fai))

    def test_one_diploid_is_one_individual(self):
        self.add("one")
        cohort = evaluate(self.data)["cohorts"][0]
        self.assertFalse(cohort["ready"])
        self.assertEqual(cohort["individuals"], ["one"])

    def test_same_species_phased_pair_ready_but_disabled(self):
        self.add("one")
        self.add("two", phased=False)
        result = evaluate(self.data)
        self.assertEqual(result["cohorts"][0]["status"], "disabled_but_eligible")
        self.assertEqual(len(result["cohorts"][0]["members"]), 4)
        self.assertTrue(all(r["qc"] == "unassessed" for r in result["assemblies"]))

    def test_fragmented_phased_member_does_not_vote(self):
        self.add("one")
        self.add("two", phased=False, fragmented=True)
        result = evaluate(self.data)
        self.assertTrue(result["cohorts"][0]["ready"])
        rows = [r for r in result["assemblies"] if r["sample"] == "two"]
        self.assertTrue(all(r["graph_member"] and not r["voter"] for r in rows))

    def test_collapsed_shortread_excluded(self):
        self.add("one", assembler="spades")
        self.add("two", assembler="spades")
        result = evaluate(self.data)
        self.assertFalse(result["cohorts"][0]["ready"])
        self.assertTrue(all(not r["representation_eligible"] for r in result["assemblies"]))

    def test_declared_haploid_shortread_allowed(self):
        self.add("one", assembler="spades", ploidy=1)
        self.add("two", assembler="spades", ploidy=1)
        self.assertTrue(evaluate(self.data)["cohorts"][0]["ready"])

    def test_missing_member_withholds_only_its_species(self):
        self.add("one", missing=("hap2",))
        self.add("two")
        self.add("three", taxid="2")
        self.add("four", taxid="2")
        cohorts = {c["taxid"]: c for c in evaluate(self.data)["cohorts"]}
        self.assertFalse(cohorts["1"]["ready"])
        self.assertTrue(cohorts["2"]["ready"])

    def test_species_do_not_pool_individuals(self):
        self.add("one")
        self.add("two", taxid="2")
        self.assertTrue(all(not c["ready"] for c in evaluate(self.data)["cohorts"]))

    def test_names_do_not_collide_and_numeric_taxid_reference_matches(self):
        self.add("a-b")
        self.add("a_b")
        self.data["harmonization_references"] = [[1, "a_b_hap2"]]
        cohort = evaluate(self.data)["cohorts"][0]
        self.assertEqual(cohort["reference_id"], "a_b_hap2")
        names = [m["graph_name"] for m in cohort["members"]]
        self.assertEqual(len(names), len(set(names)))

    def test_duplicate_final_identity_rejected(self):
        self.add("one")
        self.data["observed"].append(self.data["observed"][0])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            evaluate(self.data)


if __name__ == "__main__":
    unittest.main()