import sys
import tempfile
import unittest
from collections import defaultdict, Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_controls import sam_measure, select_controls, scan_contacts, compare_controls, peer_relationship, continuity_grid, nearby_hypotheses


class Controls(unittest.TestCase):
    def test_molecule_coverage_not_record_count(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'sam'
            line='read\t0\tscaffold\t1\t60\t4000M\t*\t0\t0\t*\t*\tNM:i:0\n'
            path.write_text(line+line+'deletion\t0\tscaffold\t1\t60\t1500M1000D1500M\t*\t0\t0\t*\t*\tNM:i:0\n')
            value=sam_measure(path,1500,2500)
            self.assertEqual(value['spanning'],1)
            self.assertEqual(value['left_molecules'],2)
            self.assertEqual(value['right_molecules'],2)
            self.assertEqual(value['left_covered_fraction'],0)

    def test_controls_exclude_other_suspect_boundaries_and_overlap(self):
        intervals={'a':dict(scaffold='s',lo=1000000,hi=1000100),
                   'b':dict(scaffold='s',lo=2000000,hi=2000100)}
        gaps=[dict(scaffold='s',lo=x,hi=x+100) for x in (1000000,1100000,2000000,3000000,3100000,4000000)]
        controls,links=select_controls(intervals,gaps,{'s':5000000})
        self.assertEqual({r['lo'] for r in controls.values()},{3000000,4000000})
        self.assertEqual(len(links['a']),2)

    def test_per_library_lift_is_one_based_and_does_not_pool(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'pairs'
            path.write_text('A_read\tsource\t8\tsource\t13\t+\t-\tUU\nB_read\tsource\t1\tsource\t9\t+\t-\tUU\n')
            placements={'source':[('s',0,20,0,20,'-')]}
            regions={'gap':dict(scaffold='s',lo=8,hi=12)}
            counts,totals,audit=scan_contacts(path,placements,regions,[dict(qname_prefix='A_',library_id='A'),dict(qname_prefix='B_',library_id='B')],flank=8)
            self.assertEqual(counts['A','gap']['cross'],1)
            self.assertEqual(counts['B','gap']['cross'],0)
            self.assertEqual(dict(totals),{'A':1,'B':1})

    def test_both_libraries_and_matched_controls_required(self):
        hifi=dict(left_molecules=15,right_molecules=15,left_covered_fraction=1,right_covered_fraction=1,spanning=0)
        focal=dict(key='gap',hifi=hifi)
        controls=[dict(key='c'+str(i),hifi=dict(hifi,spanning=3),peer_continuous_individuals=3) for i in range(5)]
        counts=defaultdict(Counter)
        for lib in ('A','B'):
            counts[lib,'gap'].update(left_within=1000,right_within=1000,cross=0)
            for control in controls:counts[lib,control['key']].update(left_within=1000,right_within=1000,cross=100)
        value=compare_controls(focal,controls,counts,{'A','B'})
        self.assertTrue(value['hic_support_loss'])
        counts['B','gap']['left_within']=0
        self.assertFalse(compare_controls(focal,controls,counts,{'A','B'})['hic_support_loss'])
        controls[0]['hifi']['spanning']=0
        self.assertFalse(compare_controls(focal,controls,counts,{'A'})['matched_controls_pass'])

    def test_gapped_paf_without_aligned_anchor_bases_is_uninformative(self):
        hit=['q','200100','0','200100','+','chr','10000000','0','100000','100000','200100','60','cg:Z:50000M100100I50000M']
        self.assertEqual(peer_relationship([hit],100000,100100)[0],'uninformative')

    def test_nearby_gap_is_discovered_but_not_silently_snapped(self):
        intervals={'transition':dict(scaffold='s',lo=1000000,hi=1010000)}
        gaps=[dict(scaffold='s',lo=1110000,hi=1110100)]
        result=nearby_hypotheses(intervals,gaps)
        self.assertEqual(len(result),1)
        self.assertEqual(intervals['transition']['hi'],1010000)
        self.assertEqual(next(iter(result.values()))['role'],'review_hypothesis')

    def test_grid_uses_narrow_spans_instead_of_broad_interval_bridges(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'sam'
            rows=[]
            for start in (0,2000,4000):
                for n in range(2):rows.append('r%d_%d\t0\ts\t%d\t60\t5000M\t*\t0\t0\t*\t*\tNM:i:0\n'%(start,n,start+1))
            path.write_text(''.join(rows))
            grid=continuity_grid(path,1500,7500)
            self.assertEqual(grid['minimum_molecules'],2)
            self.assertEqual(sam_measure(path,1500,7500)['spanning'],0)


if __name__=='__main__':unittest.main()
