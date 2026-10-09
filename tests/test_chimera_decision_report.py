import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_decisions import advise, assessment_label, decision_card, peer_counts


def track(sample,left='chr9',right='chr7',qualified=True):
    return dict(peer=sample+'_hap1',sample=sample,auto_evidence=True,
                relationship=('same_chromosome' if left==right else 'different_chromosomes') if qualified else 'uninformative',
                left=dict(chrom=left),right=dict(chrom=right))


class DecisionReport(unittest.TestCase):
    def fixture(self):
        row=dict(id='C02',assembly='assembly_hap1',scaffold='s',source_candidate_id='gap',assessment_status='review_required',
                 detected_transition='chr9 → chr7',chromosome_context='chr9 → chr7',cut_bp='325',gap_start='225',gap_end='325',
                 review_start='225',review_end='325',report_path='assembly_hap1.review/report.md#candidate-c02')
        measurement=dict(verified_gap=True,native_continuity=False,chromosome_tracks=[track('A'),track('B')])
        return row,measurement

    def test_gap_advice_requires_matching_chromosome_pair_not_any_different_pair(self):
        row,m=self.fixture()
        self.assertEqual(advise(row,m,'self')[0]['suggested_action'],'CONSIDER_GAP_CUT')
        m['chromosome_tracks'][1]=track('B',right='chr3')
        self.assertEqual(advise(row,m,'self')[0]['suggested_action'],'KEEP_CONFLICTING_EVIDENCE')

    def test_two_haplotypes_of_one_individual_and_sister_do_not_satisfy_two_individuals(self):
        row,m=self.fixture();m['chromosome_tracks']=[track('A'),track('A'),track('self')]
        self.assertEqual(peer_counts(row,m,'self')['matches'],['A'])
        self.assertEqual(advise(row,m,'self')[0]['suggested_action'],'KEEP_INSUFFICIENT_EVIDENCE')

    def test_molecule_continuity_changes_cut_advice(self):
        row,m=self.fixture();m.update(hifi_informative=True,hifi_spanning_molecules=4)
        self.assertEqual(advise(row,m,'self')[0]['suggested_action'],'KEEP_CONFLICTING_EVIDENCE')

    def test_transition_range_is_never_presented_as_exact_cut(self):
        row,m=self.fixture();row['cut_bp']='';m['verified_gap']=False
        text='\n'.join(decision_card([row],{'gap':m},{'sample':'self'}))
        self.assertIn('Keep for now — locate the break first',text)
        self.assertIn('225–325 bp; exact cut unknown',text)
        self.assertNotIn('**Cut option:**',text)

    def test_gap_card_shows_support_opposition_limit_and_exact_row(self):
        row,m=self.fixture();m['chromosome_tracks'].append(track('C',left='chr9',right='chr9',qualified=False))
        text='\n'.join(decision_card([row],{'gap':m},{'sample':'self'}))
        self.assertIn('use row C02 at **325 bp**',text)
        self.assertIn('**Evidence for cutting**',text)
        self.assertIn('**Evidence for keeping**',text)
        self.assertIn('Lower-confidence alignments place both sides on one chromosome in C',text)
        self.assertNotIn('selected=YES',text)

    def test_rejected_signal_is_keep_and_incomplete_assembly_is_not_clean(self):
        row,m=self.fixture();row['assessment_status']='retention_supported'
        self.assertEqual(advise(row,m,'self')[0]['suggested_action'],'KEEP_REJECTED_SIGNAL')
        self.assertEqual(assessment_label({'assessment_status':'assessed'},[]),'No potential chimeras detected')
        self.assertEqual(assessment_label({'assessment_status':'chromosome scope unresolved'},[]),'No chromosome-scale scaffolds identified')
        self.assertIn('Not assessed',assessment_label({'assessment_status':'unavailable'},[]))


if __name__=='__main__':unittest.main()
