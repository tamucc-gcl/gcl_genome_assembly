import unittest,sys,copy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_markdown import assess_transitions

class TransitionAssessment(unittest.TestCase):
    def fixture(self,chroms):
        rows=[dict(id='C01',source_candidate_id='a',scaffold='s',review_start=0,review_end=100,chromosome_context='Unresolved from qualified local alignments',review_priority='Investigate'),dict(id='C02',source_candidate_id='b',scaffold='s',review_start=400,review_end=410,chromosome_context='chr9 → chr7',review_priority='Prioritize gap-cut review')]
        tiles=[dict(lo=100+i*100,hi=200+i*100,chrom=c) for i,c in enumerate(chroms)]
        ms={'a':dict(packet_interval_id='p',chromosome_tracks=[dict(peer='P',sample='peer',auto_evidence=True,bins=tiles)]),'b':dict(verified_gap=True,source_candidates=['p'])}
        return rows,ms,dict(sample='self',intervals={'p':dict(start=0)}),{'a':dict(left_chrom='chr9',right_chrom='chr7'),'b':{}}
    def test_localized_transition_groups_only_with_measured_bridge(self):
        args=self.fixture(['chr9','chr7','chr7']);assess_transitions(*args)
        self.assertEqual(args[0][0]['assessment_status'],'supporting_measurement')
        self.assertEqual(args[0][0]['preferred_candidate'],'C02')
        self.assertEqual(args[0][0]['chromosome_context'],'chr9 → chr7')
        self.assertIn('transition_bridge_assessment',args[1]['b'])
    def test_reversal_prevents_precedence(self):
        args=self.fixture(['chr7','chr9','chr7'])
        for r in args[0]:
            r['review_start']*=100;r['review_end']*=100
        for b in args[1]['a']['chromosome_tracks'][0]['bins']:
            b['lo']*=100;b['hi']*=100
        assess_transitions(*args)
        self.assertEqual(args[0][0]['assessment_status'],'review_required')
        self.assertEqual(args[0][1]['bridge_status'],'reversal_or_other_chromosome_observed')
    def test_missing_bridge_does_not_merge(self):
        args=self.fixture([]);assess_transitions(*args)
        self.assertEqual(args[0][0]['assessment_status'],'review_required')
        self.assertEqual(args[0][1]['bridge_status'],'insufficient_bridge_observability')
    def test_retention_is_separate_assessment(self):
        args=self.fixture([]);args[0][1]['review_priority']='Evidence favors retaining this sampled boundary'
        assess_transitions(*args)
        self.assertEqual(args[0][1]['assessment_status'],'retention_supported')
