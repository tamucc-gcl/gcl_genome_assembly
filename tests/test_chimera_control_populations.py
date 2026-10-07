import sys,unittest
from pathlib import Path
from collections import defaultdict,Counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_controls import continuous_controls,compare_controls,peer_anchor_trials,anchor_rejection
class Populations(unittest.TestCase):
    def test_continuous_sites_exclude_gap_and_candidate_flanks(self):
        controls,links=continuous_controls({'x':dict(scaffold='s',lo=1000000,hi=1000100)},[dict(scaffold='s',lo=2000000,hi=2000100)],{'s':4000000})
        self.assertTrue(controls)
        for c in controls.values():
            self.assertEqual(c['lo'],c['hi']);self.assertEqual(c['role'],'continuous_control')
            self.assertGreater(abs(c['lo']-1000000),249999)
            self.assertGreater(abs(c['lo']-2000000),249999)
    def test_gap_controls_do_not_require_bridges_but_require_peer_continuity(self):
        h=dict(left_molecules=15,right_molecules=15,left_covered_fraction=1,right_covered_fraction=1,spanning=0)
        focal=dict(key='f',hifi=h)
        controls=[dict(key='s%d'%i,role='continuous_control',hifi=dict(h,spanning=2)) for i in range(5)]
        controls += [dict(key='g%d'%i,role='gap_control',hifi=dict(h,spanning=0),peer_continuous_individuals=3) for i in range(5)]
        counts=defaultdict(Counter)
        counts['A','f'].update(left_within=1000,right_within=1000,cross=0)
        for c in controls:counts['A',c['key']].update(left_within=1000,right_within=1000,cross=100)
        self.assertTrue(compare_controls(focal,controls,counts,{'A'})['hic_support_loss'])
        controls[-1]['peer_continuous_individuals']=0
        self.assertFalse(compare_controls(focal,controls,counts,{'A'})['hic_support_loss'])
    def test_small_unique_anchors_recover_with_explicit_failures(self):
        hits=[['q','200000','90000','110000','+','chr','10000000','90000','110000','20000','20000','60','cg:Z:20000M']]
        selected,trials=peer_anchor_trials(hits,100000,100000,200000)
        self.assertEqual(selected[0],'continuous_context')
        self.assertTrue(any(t['left_reason']=='insufficient_anchor_coverage' for t in trials))
        self.assertEqual(anchor_rejection([],0,5000),'no_alignment')
    def test_competing_targets_are_rejected(self):
        a=['q','200000','0','200000','+','one','10000000','0','200000','200000','200000','60','cg:Z:200000M']
        b=list(a);b[5]='two';b[11]='0'
        self.assertEqual(peer_anchor_trials([a,b],100000,100000,200000)[0][0],'uninformative')
