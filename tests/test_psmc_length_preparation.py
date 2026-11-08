import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'scripts/comparisons'))
spec=importlib.util.spec_from_file_location('length_prepare',REPO/'scripts/comparisons/prepare_psmc_length_tests.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


class LengthPreparation(unittest.TestCase):
    def test_only_threshold_and_output_location_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'full';control=source/'asm20'
            (control/'results/post_assembly/psmc/1/individual/genotypes').mkdir(parents=True)
            (source/'comparison.json').write_text(json.dumps(dict(sample='individual',source='/mask/source')))
            (control/'exit_status.txt').write_text('0')
            retained=root/'retained';retained.write_text('retained')
            params=dict(psmc_assembly_preset='asm20',psmc_assembly_min_alignment_bp=50000,
                psmc_final_eligibility=str(retained),psmc_manifest=str(retained),psmc_bootstraps=10,
                psmc_pattern='4+25*2+4+6',psmc_indel_padding=5,outdir='/original')
            (control/'params.json').write_text(json.dumps(params))
            job=dict(sample='individual',mode='assemblies',taxid='1',assemblies=[str(retained)]*2,reference=dict(sha256='unchanged'))
            (control/'results/post_assembly/psmc/psmc_plan.json').write_text(json.dumps(dict(jobs=[job])))
            (control/'results/post_assembly/psmc/1/individual/genotypes/dipcall_settings.json').write_text('{"assembly_preset":"asm20"}')
            with patch.object(mod,'hifi_bam',return_value=retained):specification=mod.prepare(source,root/'out',root/'hifi')
            for row in specification['trials']:
                actual=json.loads((root/'out'/row['label']/'params.json').read_text())
                expected=dict(params,psmc_assembly_min_alignment_bp=row['min_alignment_bp'],outdir=str(root/'out'/row['label']/'results'))
                self.assertEqual(actual,expected)
            self.assertEqual(specification['run_directories']['asm20_50kb'],str(control))


if __name__=='__main__':unittest.main()
