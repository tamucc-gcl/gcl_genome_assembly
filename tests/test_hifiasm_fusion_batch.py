"""Verify treatment isolation and avoid double-counting alignment span."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/comparisons'))
from run_hifiasm_fusion_batch import command, composition, TREATMENTS


class BatchTests(unittest.TestCase):
    def test_library_pairs_and_postjoin_are_only_deliberate_changes(self):
        libs={'A':('a1','a2'),'B':('b1','b2')}
        for label,keys,post in TREATMENTS:
            cmd=command('fish',48,'hifi',[libs[k] for k in keys],post)
            self.assertEqual(cmd[cmd.index('-u')+1],str(post))
            self.assertEqual(cmd[-1],'hifi')
            if keys:
                self.assertEqual(cmd[cmd.index('--h1')+1],','.join(libs[k][0] for k in keys))
                self.assertEqual(cmd[cmd.index('--h2')+1],','.join(libs[k][1] for k in keys))
            else:
                self.assertNotIn('--h1',cmd)
            self.assertNotIn('--dual-scaf',cmd)

    def test_composition_union_and_ambiguity(self):
        h=dict(query='c',target='r',query_length=10000,query_start=0,query_end=6000,mapq=60,strand='+')
        rows=composition([h,dict(h,query_start=4000,query_end=8000),dict(h,target='repeat')])
        self.assertEqual(rows[0]['union_query_span_bp'],8000)
        self.assertEqual(len(rows),2)  # competing references remain explicit


if __name__=='__main__': unittest.main()
