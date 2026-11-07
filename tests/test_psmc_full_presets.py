"""Checks that full comparisons preserve audited scope and select the real route."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'scripts/comparisons'))
spec=importlib.util.spec_from_file_location('prepare_psmc_full_presets',REPO/'scripts/comparisons/prepare_psmc_full_presets.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


class FullPresetPreparation(unittest.TestCase):
    def fixture(self,root):
        baseline=root/'baseline';(baseline/'provenance').mkdir(parents=True)
        (baseline/'provenance/run_provenance.json').write_text(json.dumps(dict(context=dict(parameters=dict(
            psmc_iterations=25,psmc_pattern='4+25*2+4+6',psmc_mode='hifi',psmc_mutation_rate=1e-8,
            psmc_bootstraps=100,psmc_min_contig_bp=999,unsupported_old_flag=True)))))
        rows=[]
        for h in ['hap1','hap2']:
            fa=root/(h+'.fa');fa.write_text('>chr1\nA\n');fai=root/(h+'.fai');fai.write_text('chr1\t1\n')
            rows.append(dict(sample='individual',haplotype=h,fasta=str(fa),fai=str(fai)))
        (baseline/'assembly_eligibility.json').write_text(json.dumps(dict(assemblies=rows)))
        source=root/'source';(source/'reports').mkdir(parents=True)
        (source/'reports/status.json').write_text('{"status":"PASS"}')
        data=dict(mode='alignment_presets',baseline=str(baseline),sample='individual',controls=[dict(
            assemblies=[r['fasta'] for r in rows],filters=dict(min_contig_bp=1000000,min_gq=20),
            include_bed='/explicit/scope.bed',exclude_bed=None)])
        (source/'manifest.json').write_text(json.dumps(data));return source

    def test_preserve_filters_fits_and_scope_but_select_consistent_unscaled_assembly_trials(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=self.fixture(root);out=root/'out'
            with patch.object(sys,'argv',['prepare','--source',str(source),'--out',str(out)]):mod.main()
            for preset in ['asm10','asm20']:
                p=json.loads((out/preset/'params.json').read_text())
                self.assertEqual(p['psmc_assembly_preset'],preset)
                self.assertEqual(p['psmc_mode'],'assemblies')
                self.assertEqual(p['psmc_min_contig_bp'],1000000)
                self.assertEqual(p['psmc_pattern'],'4+25*2+4+6')
                self.assertEqual(p['psmc_bootstraps'],10)
                self.assertIsNone(p['psmc_mutation_rate']);self.assertIsNone(p['psmc_generation_time'])
                self.assertNotIn('unsupported_old_flag',p)
            self.assertIn('/explicit/scope.bed',(out/'individual.tsv').read_text())

    def test_reject_changed_assembly_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=self.fixture(root);out=root/'out'
            path=source/'manifest.json';data=json.loads(path.read_text());data['controls'][0]['assemblies'].reverse();path.write_text(json.dumps(data))
            with patch.object(sys,'argv',['prepare','--source',str(source),'--out',str(out)]):
                with self.assertRaisesRegex(ValueError,'differ from audited'):mod.main()
            self.assertFalse(out.exists())


if __name__=='__main__':unittest.main()
