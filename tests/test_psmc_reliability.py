import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('reliability',Path(__file__).resolve().parents[1]/'scripts/comparisons/psmc_reliability.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
try:import numpy as np
except ImportError:np=None


class Reliability(unittest.TestCase):
    @unittest.skipUnless(np is not None,'Requires numpy')
    def test_scope_and_unphased_genotype_classes(self):
        a=np.frombuffer(b'AMRNCNN',dtype='S1');b=np.frombuffer(b'AMACNGN',dtype='S1')
        scope=np.array([True]*6+[False]);rows,candidates=mod.compare_block(a,b,scope)
        q={r['region']:r for r in rows}
        self.assertEqual(q['shared']['bp'],3);self.assertEqual(q['shared']['discordant'],1)
        self.assertEqual(q['shared']['het_vs_hom'],1);self.assertEqual(q['shared']['left_hets'],2)
        self.assertEqual(q['shared']['right_hets'],1);self.assertEqual(q['left_only']['bp'],1)
        self.assertEqual(q['right_only']['bp'],2);self.assertEqual(q['neither']['bp'],0)
        self.assertEqual(int(candidates['shared_discordant'].sum()),1)

    def test_read_support_requires_observability_and_both_het_alleles(self):
        self.assertEqual(mod.support('M',dict(A=4,C=4)),'insufficient_reads')
        self.assertEqual(mod.support('M',dict(A=9,C=1)),'not_supported')
        self.assertEqual(mod.support('M',dict(A=8,C=2)),'supported')
        self.assertEqual(mod.support('M',dict(A=5,C=4,G=2)),'not_supported')
        self.assertEqual(mod.support('A',dict(A=9,C=1)),'supported')
        self.assertEqual(mod.support('N',dict(A=100)),'not_called')


if __name__=='__main__':unittest.main()
