"""General junction evidence regressions, without dataset-specific locations."""
from pathlib import Path
import sys
import unittest
import tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
import misassembly_junction as j
from misassembly_enrich import attach_summaries
from misassembly_report import junction_card

class JunctionTests(unittest.TestCase):
    def test_observable_zero_uses_distance_and_depth_controls(self):
        matrix=np.full((80,80),3,dtype=np.uint32);np.fill_diagonal(matrix,100)
        matrix[:34,34:]=0;matrix[34:,:34]=0
        result=j.distance_matched_contacts(matrix,[30,31],[35,36],{'a':set(range(34)),'b':set(range(34,80))})
        self.assertEqual(result['calibration_status'],'measured')
        self.assertEqual(result['observed_expected'],0)
        self.assertAlmostEqual(result['expected_cross'],12)
        self.assertEqual(result['matched_candidate_pairs'],4)

    def test_low_identity_or_heavily_clipped_bridge_is_not_credible(self):
        row=dict(mapq=60,same_haplotype_alternatives=0,identity=.99,aligned=15000,clipped=0)
        self.assertTrue(j.credible_bridge(row))
        self.assertFalse(j.credible_bridge(dict(row,identity=.91)))
        self.assertFalse(j.credible_bridge(dict(row,clipped=12000)))
        self.assertFalse(j.credible_bridge(dict(row,same_haplotype_alternatives=1)))

    def test_missing_controls_are_not_zero_expected(self):
        matrix=np.eye(4,dtype=np.uint32)*100
        result=j.distance_matched_contacts(matrix,[0],[3],{'a':{0},'b':{3}})
        self.assertIsNone(result['expected_cross'])
        self.assertEqual(result['calibration_status'],'insufficient_matched_controls')

    def test_unobservable_anchor_is_not_calibrated(self):
        matrix=np.ones((40,40),dtype=np.uint32);np.fill_diagonal(matrix,100);matrix[10,10]=0
        result=j.distance_matched_contacts(matrix,[10],[20],{'a':set(range(40))})
        self.assertIsNone(result['observed_expected'])

    def test_depth_matching_rejects_incompatible_coverage(self):
        matrix=np.ones((40,40),dtype=np.uint32);np.fill_diagonal(matrix,100)
        matrix[10,10]=10000;matrix[20,20]=10000
        result=j.distance_matched_contacts(matrix,[10],[20],{'a':set(range(40))})
        self.assertIsNone(result['expected_cross'])

    def test_positive_arm_bridges_are_cut_counterevidence(self):
        event=dict(id='C80',lo=100,hi=200,recommendation='SUSPECT',junction_support={'unique_arm_bridges':3})
        j.refine_advice(event)
        self.assertEqual(event['recommendation'],'REVIEW');self.assertIn('unique sequence',event['reason'])

    def test_support_at_midpoint_does_not_veto_a_different_cut_site(self):
        event=dict(id='C80',lo=100,hi=200,recommendation='CUT',junction_support={'unique_arm_bridges':3},recommended_cut={'id':'G1'})
        j.refine_advice(event);self.assertEqual(event['recommendation'],'CUT')

    def test_no_read_reach_and_unknown_counts_are_explicit(self):
        event=dict(id='C80',junction_evidence=dict(sites=[dict(id='C80T',point=1000,credible_bridges=None,unique_arm_bridges=None,
            minimum_anchor_span_bp='100000',local_reads_long_enough='0',largest_independent_opposite_island_bp=10,contacts=[])]))
        text='\n'.join(junction_card(event,'test'))
        self.assertIn('no local read long enough',text);self.assertIn('| ? | ? |',text)
        self.assertNotIn('bridge_placements.tsv',text)

    def test_gap_does_not_get_false_negative_uninterrupted_read_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for branch in ('labels','hifi','hic'):(root/branch).mkdir()
            (root/'labels/opposite_side_islands.tsv').write_text('site\teligible\tsample\tlargest_unique_block_bp\n')
            (root/'hifi/unique_anchor_reach.tsv').write_text('site\tminimum_unique_anchor_span_bp\tlocal_reads_long_enough\tlocal_read_p95\nG1\t2000\t5\t15000\n')
            data=dict(sample='a',events=[dict(id='C80',junction_sites=[dict(id='G1',point=1000,kind='gap_option')],cut_options=[dict(id='G1')])])
            attach_summaries(data,root)
            site=data['events'][0]['junction_evidence']['sites'][0]
            self.assertIsNone(site['credible_bridges']);self.assertIsNone(site['unique_arm_bridges'])

    def test_windows_include_every_option(self):
        event=dict(id='C80',lo=100,hi=200,alternatives=[dict(id='C80R01',lo=300,hi=400,probe_bp=350)],cut_options=[])
        sites,lo,hi=j.windows(event,1000,padding=50)
        self.assertEqual([s['point'] for s in sites],[150,350]);self.assertEqual((lo,hi),(50,450))

if __name__=='__main__':unittest.main()
