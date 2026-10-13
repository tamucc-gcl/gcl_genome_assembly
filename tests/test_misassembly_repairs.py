import csv
import json
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'py_scripts'))
from misassembly_core import discover,gap_association,flank_support,support_grid,contact_result,matched_control_sites,recommend


class EvidenceRepairTests(unittest.TestCase):
    def test_missing_consensus_bins_do_not_keep_opposed_single_peer_fusion(self):
        def blocks(left):return [dict(lo=0,hi=2000000,chrom=left),dict(lo=2200000,hi=4200000,chrom='chr5')]
        peers=[(dict(id=p,sample=p,eligible=True),{'s':blocks('chr2' if p=='a' else 'chr5')}) for p in ('a','b','c')]
        events,audit=discover(peers,'focal')
        self.assertEqual(events,[])
        rejection=[r for r in audit if r.get('kind')=='rejected_peer_transition'][0]
        self.assertEqual(rejection['peers']['continuous'],['b','c'])
        self.assertEqual(rejection['peers']['separate'],['a'])
        # Independent disagreement is retained for review, not majority-suppressed.
        peers.append((dict(id='d',sample='d',eligible=True),{'s':blocks('chr2')}))
        self.assertEqual(len(discover(peers,'focal')[0]),1)

    def test_failed_fine_refinement_keeps_coarse_resolution(self):
        blocks=[dict(lo=0,hi=2000000,chrom='chr1'),dict(lo=2100000,hi=4100000,chrom='chr2')]
        peers=[(dict(id=p,sample=p,eligible=True),{'s':blocks}) for p in ('a','b')]
        coarse=[dict(b,aligned_bp=b['hi']-b['lo']) for b in blocks]
        with patch('misassembly_core.consensus_regions',side_effect=[coarse,[]]):
            events,_=discover(peers,'focal')
        self.assertEqual(events[0]['localization_resolution_bp'],100000)

    def event(self,chrom=None):
        tracks=[dict(peer=p,sample=p,eligible=True,blocks=[dict(lo=1000,hi=3000,chrom=chrom)] if chrom else []) for p in ('a','b')]
        return dict(lo=1000,hi=2000,left='chr1',right='chr2',peer_tracks=tracks)

    def test_association_missing_is_not_contradiction(self):
        event=self.event()
        result=gap_association(event,dict(lo=2500,hi=2600),'focal')
        self.assertEqual(result['status'],'unresolved')
        self.assertEqual(result['opposing_individuals'],[])
        self.assertEqual(gap_association(self.event('chr3'),dict(lo=2500,hi=2600),'focal')['status'],'contradicted')
        self.assertEqual(gap_association(self.event('chr2'),dict(lo=2500,hi=2600),'focal')['status'],'compatible')

    def test_association_respects_localization_resolution(self):
        event=dict(self.event(),localization='Independent-individual regional chromosome transition',localization_resolution_bp=1000)
        self.assertEqual(gap_association(event,dict(lo=2225,hi=2325),'focal')['status'],'confirmed')
        self.assertEqual(gap_association(event,dict(lo=5000,hi=5100),'focal')['status'],'unresolved')

    def test_adjacent_neighborhood_can_be_observed_despite_sparse_endpoint(self):
        with tempfile.TemporaryDirectory() as temp:
            sam=Path(temp)/'reads.sam'
            rows=[]
            for i in range(8):
                rows += [f'l{i}\t0\ts\t1\t60\t19900M\t*\t0\t0\t*\t*\tNM:i:0\n',f'r{i}\t0\ts\t20101\t60\t20000M\t*\t0\t0\t*\t*\tNM:i:0\n']
            sam.write_text(''.join(rows))
            flanks=flank_support(sam,20000,20100,50000)
            self.assertTrue(all(f['observable'] for f in flanks))
            self.assertEqual(flanks[0]['minimum'],0)
            self.assertEqual(flanks[0]['median'],8)
            self.assertEqual(flanks[0]['supported_fraction'],.9)
            sparse=flank_support(sam,500,600,50000)
            self.assertFalse(sparse[0]['observable'])

    def test_last_unsupported_probe_is_a_position_not_zero_length_interval(self):
        with tempfile.TemporaryDirectory() as temp:
            sam=Path(temp)/'reads.sam'
            sam.write_text(''.join(f'r{i}\t0\ts\t1\t60\t3999M\t*\t0\t0\t*\t*\tNM:i:0\n' for i in range(2)))
            grid=support_grid(sam,1500,3000)
            self.assertEqual(grid['weak_runs'],[])
            self.assertEqual(grid['weak_points'],[dict(position=3000,molecules=0,boundary_censored=True)])
            self.assertEqual(grid['minimum'],0)

    def test_control_audit_retains_context_depth_and_observability_exclusions(self):
        focal=dict(cross=0,left_within=157,right_within=4819)
        controls=[dict(id='context',cross=100,left_within=500,right_within=5000,context_rejection='insufficient_independent_same_chromosome_context'),
                  dict(id='depth',cross=100,left_within=20000,right_within=5000),dict(id='low',cross=0,left_within=0,right_within=5000),
                  dict(id='match',cross=100,left_within=500,right_within=5000)]
        result=contact_result(focal,controls)
        self.assertEqual(result['tested_controls'],4)
        self.assertEqual(result['control_ids'],['match'])
        self.assertEqual([r['rejection'] for r in result['control_audit']],['insufficient_independent_same_chromosome_context','flank_depth_mismatch','insufficient_within_side_contacts','accepted'])
        self.assertFalse(result['calibrated'])
        self.assertEqual(result['calibration_failure'],['fewer_than_5_matched_controls'])
        self.assertEqual(result['control_audit'][0]['left_within'],500)

    def test_control_pool_represents_multiple_scaffolds(self):
        gaps=[dict(scaffold=sc,lo=x,hi=x+100) for sc in ('a','b') for x in (1000000,3000000,5000000,7000000)]
        controls,links=matched_control_sites({'G':dict(scaffold='a',lo=2000000,hi=2000100,cut_bp=2000100)},gaps,{'a':10000000,'b':10000000},count=2)
        self.assertEqual({c['scaffold'] for c in controls.values()},{'a','b'})
        self.assertEqual(len(links['G']),2)

    def test_two_library_loss_and_read_continuity_is_explicit_conflict(self):
        event=dict(peers=dict(separate=['a','b'],continuous=[],conflicting=[]),hifi=dict(state='continuous'),hifi_strict=dict(state='continuous'),contacts=[dict(library='A',calibrated=True,loss=True),dict(library='B',calibrated=True,loss=True)])
        self.assertEqual(recommend(event)[0],'REVIEW')
        self.assertEqual(len(event['evidence_conflicts']),2)
        event['contacts'].pop()
        self.assertEqual(recommend(event)[0],'REVIEW')
        self.assertEqual(len(event['evidence_conflicts']),1)

    def test_one_library_can_support_a_cut_without_auto_selection(self):
        peers=dict(separate=['a','b'],continuous=[],conflicting=[])
        option=dict(id='G1',peers=peers,spanning=0,spanning_all_primary=0,flank_support=[dict(observable=True),dict(observable=True)],association=dict(status='confirmed'),contacts=[dict(library='A',calibrated=True,loss=True)])
        event=dict(peers=peers,hifi=dict(state='patchy'),cut_options=[option])
        self.assertEqual(recommend(event)[0],'CUT')
        self.assertEqual(event['recommended_cut']['id'],'G1')
        option['contacts'].append(dict(library='B',calibrated=False,loss=False))
        self.assertEqual(recommend(event)[0],'CUT')
        option['contacts'][-1]['calibrated']=True
        self.assertNotEqual(recommend(event)[0],'CUT')
        self.assertTrue(any('disagreement' in v for v in option['unmet_cut_prerequisites']))
        option['contacts']=[]
        self.assertNotEqual(recommend(event)[0],'CUT')
        self.assertTrue(any('No available Hi-C' in v for v in option['unmet_cut_prerequisites']))

    def test_stricter_read_check_can_prevent_retention(self):
        event=dict(peers=dict(separate=['a','b'],continuous=[],conflicting=[]),hifi=dict(state='continuous'),hifi_strict=dict(state='patchy'))
        self.assertEqual(recommend(event)[0],'REVIEW')
        self.assertIn('MAPQ-30',event['reason'] if 'reason' in event else event['evidence_conflicts'][0])
        event['hifi_strict']=dict(state='continuous')
        self.assertEqual(recommend(event)[0],'RETAIN')

    def test_retained_transitions_stay_visible_in_overview_and_reports(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);packet=root/'packet';packet.mkdir()
            event=dict(id='C01',scaffold='s',lo=100,hi=200,left='chr1',right='chr2',recommendation='RETAIN',reason='Local reads support continuity.',peers=dict(separate=['a','b'],continuous=[],conflicting=[]),cut_options=[],peer_tests=[])
            data=dict(assembly='fixture',sample='fixture',assessment_sha256='fixture',status='assessed',mapping=dict(status='unavailable'),events=[event])
            (packet/'evidence.json').write_text(json.dumps(data))
            subprocess.run([sys.executable,str(ROOT/'py_scripts/misassembly_report.py'),'--packet',str(packet)],cwd=root,check=True,capture_output=True)
            main=(root/'README.md').read_text(encoding='utf-8');detail=(root/'assemblies/fixture/report.md').read_text(encoding='utf-8')
            self.assertIn('1 chromosome transition(s) detected; 0 need review; 1 provisionally keep intact',main)
            self.assertIn('fixture C01',main)
            self.assertIn('### C01',detail)
            self.assertNotIn('Boundaries retained after assessment',detail)

    def test_empty_packet_publishes_only_existing_evidence_links(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);packet=root/'packet';packet.mkdir()
            data=dict(assembly='empty',sample='empty',assessment_sha256='fixture',status='assessed',mapping=dict(status='not_needed'),events=[])
            (packet/'evidence.json').write_text(json.dumps(data))
            subprocess.run([sys.executable,str(ROOT/'py_scripts/misassembly_report.py'),'--packet',str(packet)],cwd=root,check=True,capture_output=True)
            report=root/'assemblies/empty/report.md'
            for path in re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf-8')):
                self.assertTrue((report.parent/path.split('#')[0]).exists(),path)
            self.assertNotIn('read_support.tsv',report.read_text(encoding='utf-8'))


if __name__=='__main__':unittest.main()
