import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from misassembly_core import regional_hypotheses,assess_alternative,rank_alternatives,recommend,clip_block,discover
from misassembly_assess import bridge_records,chromosome_colors
from misassembly_report import decisions


class RegionalEvidenceTests(unittest.TestCase):
    def event(self):
        probes=[dict(position=p,molecules=8) for p in range(0,100001,1000)]
        return dict(id='C01',scaffold='s',left='chr1',right='chr2',lo=40000,hi=41000,start=0,end=100000,
                    peer_tracks=[],peers=dict(separate=['a','b'],continuous=[],conflicting=[]),cut_options=[],
                    hifi=dict(state='continuous'),hifi_strict=dict(state='continuous'),
                    **{key:dict(probes=[dict(p) for p in probes]) for key in ('hifi_profile','hifi_strict_profile','hifi_all_primary_profile')})

    def test_adjacent_zero_region_is_assessed_and_prevents_narrow_retention(self):
        event=self.event()
        for key in ('hifi_profile','hifi_strict_profile'):
            for p in event[key]['probes']:
                if 60000<=p['position']<70000:p['molecules']=0
        event['alternatives']=regional_hypotheses(event)
        self.assertEqual(len(event['alternatives']),1)
        target=event['alternatives'][0]
        self.assertEqual((target['lo'],target['hi']),(60000,70000))
        assess_alternative(event,target,'focal')
        rank_alternatives(event)
        self.assertTrue(target['ambiguity_restores_support'])
        self.assertEqual(recommend(event)[0],'REVIEW')
        self.assertIsNone(target['cut_bp'])
        data=dict(assembly='focal_hap1',assessment_sha256='abc',events=[dict(event,recommendation='REVIEW',reason='Review measured alternatives')])
        row=decisions(data)[1]
        self.assertEqual(row['selected'],'NO')
        self.assertIsNone(row['cut_bp'])
        self.assertEqual(row['assay_probe_bp'],60000)

    def test_strict_only_sustained_loss_is_distinguished_from_missing_reads(self):
        event=self.event()
        for p in event['hifi_strict_profile']['probes']:
            if 60000<=p['position']<66000:p['molecules']=0
        target=regional_hypotheses(event)[0]
        assess_alternative(event,target,'focal')
        self.assertEqual(target['hifi']['minimum'],8)
        self.assertEqual(target['hifi_strict']['minimum'],0)
        self.assertTrue(target['ambiguity_restores_support'])

    def test_missing_profiles_do_not_fabricate_zero_support(self):
        event=self.event();event['hifi_profile']={};event['hifi_strict_profile']={}
        self.assertEqual(regional_hypotheses(event),[])

    def test_isolated_strict_dip_is_not_a_sustained_alternative(self):
        event=self.event();event['hifi_strict_profile']['probes'][60]['molecules']=0
        self.assertEqual(regional_hypotheses(event),[])

    def test_regions_are_clipped_and_limit_is_audited(self):
        event=self.event()
        for p in event['hifi_profile']['probes']:
            if p['position']%2000==0:p['molecules']=0
        candidates=regional_hypotheses(event,maximum=2)
        self.assertEqual(len(candidates),2)
        self.assertGreater(event['unassessed_alternative_count'],0)
        self.assertTrue(all(t['hi']<=event['end'] for t in candidates))

    def test_clipping_preserves_alignment_density(self):
        block=clip_block(dict(lo=0,hi=10000,chrom='chr1',aligned_bp=5000),2000,4000)
        self.assertEqual(block['aligned_bp'],1000)

    def test_actual_bridges_deduplicate_and_audit_competing_placements(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'reads.sam'
            path.write_text('r1\t0\ts\t1\t10\t10000M\t*\t0\t0\t*\t*\tNM:i:0\n'*2+
                            'r1\t256\ts\t1\t0\t10000M\t*\t0\t0\t*\t*\tNM:i:0\n'+
                            'r2\t0\ts\t1\t60\t4000M100D5900M\t*\t0\t0\t*\t*\tNM:i:100\n')
            rows=bridge_records(path,[dict(id='R01',probe_bp=4000)])
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['mapq'],10)
            self.assertTrue(rows[0]['competing_placement_in_window'])

    def test_focal_colors_are_distinct_and_stable_across_haplotypes(self):
        colors=chromosome_colors(['chr7','chr12','chr4'],['chr7','chr12'])
        reverse=chromosome_colors(['chr7','chr12'],['chr12','chr7'])
        self.assertNotEqual(colors['chr7'],colors['chr12'])
        self.assertEqual(colors['chr7'],reverse['chr7'])
        self.assertEqual(colors['chr12'],reverse['chr12'])

    def test_unlocalized_composite_is_not_retained_from_one_supported_site(self):
        event=self.event();event['localization_warning']='Flanks insufficient'
        self.assertEqual(recommend(event)[0],'REVIEW')

    def test_fine_assignment_island_cannot_override_sustained_peer_flanks(self):
        blocks=[dict(lo=0,hi=2000000,chrom='chr1'),dict(lo=3100000,hi=5000000,chrom='chr2')]
        peers=[(dict(id=name,sample=name,eligible=True),{'s':blocks}) for name in ('a','b')]
        tiny=[dict(lo=3000000,hi=3001000,chrom='chr1'),dict(lo=3002000,hi=3200000,chrom='chr2')]
        with patch('misassembly_core.consensus_regions',side_effect=[blocks,tiny]):
            events,_=discover(peers,'focal')
        self.assertEqual(len(events),1)
        self.assertEqual((events[0]['lo'],events[0]['hi']),(2000000,3100000))
        self.assertIn('localization_warning',events[0])


if __name__=='__main__':unittest.main()
