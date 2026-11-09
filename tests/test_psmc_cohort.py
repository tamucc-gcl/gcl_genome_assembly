import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'scripts/comparisons'))
spec=importlib.util.spec_from_file_location('psmc_cohort',REPO/'scripts/comparisons/psmc_cohort.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


class CohortSelection(unittest.TestCase):
    def test_one_individual_hifi_override_is_reset_but_quality_parameters_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);base=root/'baseline';(base/'provenance').mkdir(parents=True)
            params=dict(psmc_manifest='/old/one-individual.tsv',psmc_mode='hifi',psmc_assembly_preset='asm5',
                psmc_assembly_min_alignment_bp=50000,psmc_min_gq=25,psmc_bootstraps=10,psmc_mutation_rate=1e-8,psmc_generation_time=2)
            (base/'provenance/run_provenance.json').write_text(json.dumps(dict(context=dict(parameters=params))))
            (base/'assembly_eligibility.json').write_text(json.dumps(dict(assemblies=[dict(sample='B'),dict(sample='A'),dict(sample='A')])) )
            out=root/'out';mod.prepare(base,out)
            got=json.loads((out/'baseline-params.json').read_text())
            self.assertIsNone(got['psmc_manifest']);self.assertEqual(got['psmc_mode'],'auto')
            self.assertEqual(got['psmc_assembly_preset'],'asm20');self.assertEqual(got['psmc_assembly_min_alignment_bp'],10000)
            self.assertEqual(got['psmc_min_gq'],25);self.assertEqual(got['psmc_bootstraps'],10)
            self.assertIsNone(got['psmc_mutation_rate']);self.assertIsNone(got['psmc_generation_time'])
            self.assertEqual(json.loads((out/'psmc_expected.json').read_text())['samples'],['A','B'])

    def test_environment_uses_exact_successful_input_task(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'trace.tsv').write_text('name\tstatus\texit\tworkdir\nPSMC_ANALYSIS:PSMC_INPUT (sample)\tCACHED\t0\t/exact/task\n')
            (root/'execution_environments.json').write_text(json.dumps(dict(tasks=[dict(workdir='/exact/task',conda_prefixes=['/analysis/env'])])))
            self.assertEqual(mod.environment(root,'sample'),'/analysis/env')
            with self.assertRaises(ValueError):mod.environment(root,'other')


if __name__=='__main__':unittest.main()
