import ast
import importlib.util
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/comparisons'))
path=ROOT/'scripts/comparisons/run_ctlk_boundary_batch.py'
ast.parse(path.read_text(),filename=str(path),feature_version=(3,11))
spec=importlib.util.spec_from_file_location('boundary_batch',path)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


class BoundaryTests(unittest.TestCase):
    def test_unjoin_preserves_gap_bases(self):
        left,right=mod.split_record('AACNNNGTT',3,6)
        self.assertEqual(left,'AACNNN');self.assertEqual(right,'GTT')
        self.assertEqual(left+right,'AACNNNGTT')

    def test_non_gap_refuses_sequence_cut(self):
        with self.assertRaises(ValueError):mod.split_record('AACACCGTT',3,6)

    def hit(self,strand='+'):
        return dict(ts=100,te=113,qs=0,qe=12,strand=strand,tags='cg:Z:5=2I3=3D2=')

    def test_inserted_query_coordinates_do_not_get_projected(self):
        hit=self.hit()
        self.assertIsNone(mod.project_query(hit,6))
        self.assertEqual(mod.project_query(hit,8),106)

    def test_reverse_paf_projection(self):
        hit=self.hit('-')
        self.assertEqual(mod.project_query(hit,12),100)
        self.assertEqual(mod.project_query(hit,0),113)
        self.assertIsNone(mod.project_query(hit,6))

    def test_inconsistent_paf_refuses_lift(self):
        hit=self.hit();hit['te']=114
        with self.assertRaises(ValueError):mod.query_target_blocks(hit)


if __name__=='__main__':unittest.main()
