import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from misassembly_core import recommend,site_review,regional_contact_evidence,rank_alternatives,regional_hypotheses
from misassembly_report import decisions,review_locations,event_card


class SitePolicyTests(unittest.TestCase):
    def fixture(self,point=3000000):
        broad=dict(separate=['one','two','three'],continuous=[],conflicting=[])
        local=dict(separate=['one'],continuous=[],conflicting=[])
        profiles={'library':[
            dict(position=point+d,ratio=.002 if abs(d)<=100000 else .02,cross=10,left_within=1000,right_within=1000)
            for d in range(-800000,800001,100000)]}
        gap=dict(id='boundary',lo=point-50,hi=point+50,cut_bp=point+50,peers=local,spanning=0,spanning_all_primary=0,
                 flank_support=[dict(observable=True),dict(observable=True)],association=dict(status='unresolved'),
                 contacts=[dict(library='library',calibrated=False,loss=False)])
        return dict(id='event',scaffold='anonymous',left='chr41',right='chr72',lo=point-200000,hi=point+200000,
                    independent_individuals=['one','two','three'],peers=broad,cut_options=[gap],
                    hifi=dict(state='mostly_continuous'),hifi_strict=dict(state='mostly_continuous'),
                    localization_warning='Fine localization incomplete',contact_profiles=profiles,alternatives=[])

    def test_manual_rule_is_invariant_under_coordinate_and_identifier_changes(self):
        for point in (3000000,17000000,113000000):
            e=self.fixture(point);e['id']='unseen'+str(point);e['scaffold']='unrelated_scaffold'
            self.assertEqual(recommend(e)[0],'CUT')
            self.assertEqual(e['recommended_cut']['cut_bp'],point+50)
            self.assertFalse(e['recommended_cut']['automatic_cut_eligible'])

    def test_unmeasured_low_confidence_support_is_not_zero(self):
        e=self.fixture();e['cut_options'][0].pop('spanning_all_primary')
        self.assertNotEqual(recommend(e)[0],'CUT')

    def test_supported_gap_is_never_recommended_for_cut(self):
        e=self.fixture();e['cut_options'][0]['spanning']=4
        self.assertNotEqual(recommend(e)[0],'CUT')

    def test_repeat_gap_with_low_confidence_bridges_is_not_a_cut(self):
        e=self.fixture();g=e['cut_options'][0];g['spanning_all_primary']=12
        g['association']['status']='confirmed';g['peers']['separate']=['one','two']
        g['contacts'][0].update(calibrated=True,loss=True)
        self.assertNotEqual(recommend(e)[0],'CUT')
        self.assertFalse(g['automatic_cut_eligible'])

    def test_same_chromosome_context_excludes_the_site(self):
        e=self.fixture();g=e['cut_options'][0];g['peers']=dict(separate=[],continuous=['two','three'],conflicting=[])
        self.assertEqual(site_review(e,g)[0],'EXCLUDE')
        self.assertNotEqual(recommend(e)[0],'CUT')

    def test_library_disagreement_blocks_manual_recommendation(self):
        e=self.fixture();g=e['cut_options'][0]
        g['contacts']=[dict(library='a',calibrated=True,loss=True),dict(library='b',calibrated=True,loss=False)]
        self.assertNotEqual(recommend(e)[0],'CUT')

    def test_unobserved_contact_ratios_cannot_corroborate_a_gap(self):
        e=self.fixture()
        for row in e['contact_profiles']['library']:row['right_within']=0;row['ratio']=None
        self.assertNotEqual(recommend(e)[0],'CUT')

    def test_missing_contact_shoulders_do_not_fabricate_depletion(self):
        e=self.fixture();e['contact_profiles']['library']=e['contact_profiles']['library'][7:10]
        result=regional_contact_evidence(e,e['cut_options'][0])[0]
        self.assertIsNone(result['fraction']);self.assertFalse(result['depleted'])
        self.assertNotEqual(recommend(e)[0],'CUT')

    def test_multiple_qualified_gap_options_need_a_human_choice(self):
        e=self.fixture();e['cut_options'].append(copy.deepcopy(e['cut_options'][0]));e['cut_options'][1]['id']='other'
        self.assertNotEqual(recommend(e)[0],'CUT')
        self.assertNotIn('recommended_cut',e)

    def test_single_strict_probe_loss_is_not_a_likely_break(self):
        e=self.fixture();e['cut_options']=[];e.pop('localization_warning')
        e['hifi']=dict(state='continuous',minimum=10)
        e['hifi_strict']=dict(state='mostly_continuous',minimum=1,probes=[dict(position=3000000,molecules=1)])
        self.assertNotEqual(site_review(e,e)[0],'BREAK_REGION')

    def test_region_has_no_automatic_exact_cut_coordinate(self):
        e=self.fixture();e['cut_options']=[];e.pop('localization_warning')
        e['hifi']=dict(state='continuous',minimum=10)
        e['hifi_strict']=dict(state='patchy',minimum=1,longest_weak_run_bp=7000)
        verdict,reason=recommend(e);e.update(recommendation=verdict,reason=reason)
        self.assertEqual(verdict,'SUSPECT');self.assertNotIn('recommended_cut',e)
        rows=decisions(dict(assembly='unknown',assessment_sha256='digest',events=[e]))
        self.assertEqual(rows[0]['cut_bp'],'');self.assertEqual(rows[0]['selected'],'NO')

    def test_regional_assay_and_gap_share_one_editable_option(self):
        e=self.fixture();g=e['cut_options'][0]
        e['alternatives']=[dict(id='assay',lo=g['lo']-1000,hi=g['hi']+1000,kind='regional_read_loss')]
        self.assertEqual(len(review_locations(e)),1)
        self.assertEqual(g['assay_aliases'],['assay'])

    def test_excluded_same_chromosome_trough_does_not_block_retention(self):
        e=self.fixture();e['cut_options']=[];e.pop('localization_warning')
        e['alternatives']=[dict(id='ordinary_gap',lo=3500000,hi=3510000,distance_bp=500000,kind='regional_read_loss',
            peers=dict(separate=[],continuous=['one','two'],conflicting=[]),hifi=dict(minimum=0),
            hifi_strict=dict(minimum=0),flank_support=[dict(observable=True)]*2,ambiguity_restores_support=False)]
        rank_alternatives(e)
        self.assertEqual(e['alternatives'][0]['site_verdict'],'EXCLUDE')
        self.assertFalse(e['regional_concerns'])
        self.assertEqual(recommend(e)[0],'RETAIN')

    def test_scattered_strict_losses_do_not_count_as_sustained(self):
        e=self.fixture();e['cut_options']=[];e.pop('localization_warning')
        e['hifi']=dict(state='continuous',minimum=10)
        e['hifi_strict']=dict(state='patchy',minimum=1,probes=[dict(position=3000000+i*1000,molecules=1 if i%2==0 else 10) for i in range(11)])
        self.assertNotEqual(site_review(e,e)[0],'BREAK_REGION')

    def test_low_count_contact_shoulders_do_not_create_a_baseline(self):
        e=self.fixture()
        for row in e['contact_profiles']['library']:
            if abs(row['position']-3000000)>150000:row['left_within']=1
        self.assertIsNone(regional_contact_evidence(e,e['cut_options'][0])[0]['baseline'])
        self.assertNotEqual(recommend(e)[0],'CUT')

    def test_preferred_region_is_exported_without_a_cut(self):
        e=self.fixture();e['cut_options']=[];e.update(recommendation='SUSPECT',reason='region',preferred_region=dict(lo=3100000,hi=3120000))
        rows=decisions(dict(assembly='unknown',assessment_sha256='digest',events=[e]))
        self.assertEqual((rows[0]['review_start'],rows[0]['review_end']),(3100000,3120000))
        self.assertEqual(rows[0]['cut_bp'],'')

    def test_supported_bridges_are_counterevidence_in_summary(self):
        e=self.fixture();e.update(recommendation='REVIEW',reason='review');e['cut_options'][0]['spanning']=8
        lines=event_card(e,'unknown',cohort=True)
        row=next(line for line in lines if '| boundary |' in line)
        cells=row.split('|')
        self.assertNotIn('8 gap-spanning',cells[4])
        self.assertIn('8 gap-spanning',cells[5])

    def test_wide_strict_loss_does_not_swallow_separate_confident_troughs(self):
        e=self.fixture();e['start']=2000000;e['end']=4000000
        e['hifi_profile']=dict(probes=[dict(position=p,molecules=0 if p in (3000000,3001000,3010000,3011000) else 8) for p in range(2990000,3030000,1000)])
        e['hifi_strict_profile']=dict(probes=[dict(p,molecules=0) for p in e['hifi_profile']['probes']])
        candidates=regional_hypotheses(e)
        self.assertEqual(len(candidates),2)
        self.assertEqual(sorted((t['lo'],t['hi']) for t in candidates),[(3000000,3002000),(3010000,3012000)])


if __name__=='__main__':unittest.main()
