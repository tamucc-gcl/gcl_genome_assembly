"""Adversarial evidence tests; executed in the cluster batch before tools run."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/comparisons'))
spec=importlib.util.spec_from_file_location('ctlk_adjudication',ROOT/'scripts/comparisons/run_ctlk_adjudication.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


class EvidenceTests(unittest.TestCase):
    def hit(self,target,start,matches=9900,mapq=60):
        return dict(qs=0,qe=10000,matches=matches,block=10000,mapq=mapq,target=target,ts=start,strand='+')

    def test_high_mapq_does_not_override_competing_chromosome(self):
        state,_=mod.classify_anchor([self.hit('chr4',0),self.hit('chr14',0,9800)],10000)
        self.assertEqual(state,'ambiguous')

    def test_competing_repeat_on_same_chromosome_is_ambiguous(self):
        state,_=mod.classify_anchor([self.hit('chr4',0),self.hit('chr4',1000000)],10000)
        self.assertEqual(state,'ambiguous')

    def test_short_alignment_cannot_assign_anchor(self):
        h=self.hit('chr4',0);h['qe']=1000
        self.assertEqual(mod.classify_anchor([h],10000)[0],'no_full_placement')

    def test_graph_one_hop_does_not_expand_recursively(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'g.gfa';p.write_text('L\ta\t+\tb\t-\t10M\nL\tb\t+\tc\t+\t10M\n')
            self.assertEqual(mod.graph_neighborhood(p,{'a'}),{'a','b'})

    def test_longest_sam_sequence_and_reverse_orientation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'r.sam'
            p.write_text('r\t2064\tc\t1\t60\t2H2M\t*\t0\t0\tAC\tII\n'
                         'r\t16\tc\t1\t60\t4M\t*\t0\t0\tAACC\tIIII\n')
            seqs,clipped=mod.retrieve_sequences(p)
            self.assertEqual(seqs['r'],'GGTT');self.assertIn('r',clipped)

    def test_distance_exposure_respects_gap_and_orientation(self):
        self.assertEqual(mod.distance_area((0,10),(20,30),0,10),0)
        self.assertEqual(mod.distance_area((0,10),(20,30),0,100),100)
        self.assertEqual(mod.distance_area((0,10),(0,10),0,100),50)


if __name__=='__main__':unittest.main()
