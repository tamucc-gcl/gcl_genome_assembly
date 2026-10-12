import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'py_scripts'))
from misassembly_core import cut_prerequisites,table
from misassembly_report import sister_context,decisions,event_card,gap_summary


def fixture():
    peers=dict(separate=['a','b'],continuous=[],conflicting=[])
    gap=dict(id='C01G01',lo=1000,hi=1100,cut_bp=1100,peers=peers,spanning=0,
             flank_support=[dict(side=side,observable=True,minimum=5,median=5,maximum=5,supported_fraction=1) for side in ('left','right')],
             association=dict(status='unresolved',reason='Intervening assignment is incomplete.'),
             contacts=[dict(library=library,calibrated=True,loss=False,raw_counts=dict(cross=100,left_within=1000,right_within=1000),matched_controls=5) for library in ('A','B')])
    event=dict(id='C01',scaffold='s',lo=2000,hi=2100,left='chr1',right='chr2',recommendation='REVIEW',reason='Conflicting evidence.',
               peers=peers,peer_tests=[],cut_options=[gap])
    return dict(assembly='sample_hap1',sample='sample',assessment_sha256='fixture',status='assessed',mapping=dict(status='unavailable'),events=[event])


class FinishedReviewTests(unittest.TestCase):
    def test_nonloss_in_both_libraries_is_not_library_disagreement(self):
        option=fixture()['events'][0]['cut_options'][0]
        self.assertFalse(any('disagreement' in v for v in cut_prerequisites(option)))
        option['contacts'][0]['loss']=True
        self.assertTrue(any('disagreement' in v for v in cut_prerequisites(option)))

    def test_gap_row_preserves_event_and_labels_unassociated_option(self):
        data=fixture();row=decisions(data)[0]
        self.assertEqual(row['event_id'],'C01')
        self.assertEqual(row['event_assessment'],'REVIEW')
        self.assertEqual(row['event_reason'],'Conflicting evidence.')
        self.assertIn('Review',row['recommendation'])
        self.assertIn('Not suggested',row['option_advice'])
        self.assertEqual(row['option_association'],'unresolved')
        self.assertEqual(row['selected'],'NO')
        self.assertIn('nearby only; not associated',gap_summary(data['events'][0]))

    def test_zero_gap_support_is_distinguished_from_neighboring_reads(self):
        data=fixture();event=data['events'][0];gap=event['cut_options'][0]
        gap['association']=dict(status='confirmed',reason='Within measured resolution.')
        event['hifi']=dict(probes=[dict(position=2000,molecules=3)],supported_fraction=.875,minimum=0,state='patchy')
        gap['contacts']=[];gap.pop('flank_support')
        text='\n'.join(event_card(event,data['assembly']))
        self.assertIn('Neighboring-sequence support',text)
        self.assertIn('zero HiFi molecules span the actual gap',text)
        self.assertIn('1,100 bp (available; not recommended)',gap_summary(event))

    def test_counterpart_detection_and_new_manual_columns_are_rendered(self):
        first=fixture();second=fixture();second['assembly']='sample_hap2'
        second['events'][0].update(left='chr2',right='chr1',scaffold='other')
        stranger=fixture();stranger.update(assembly='unrelated',sample='another')
        self.assertEqual(len(sister_context(first['events'][0],first,[first,second,stranger])),1)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);args=[]
            for i,data in enumerate([first,second]):
                packet=root/str(i);packet.mkdir();(packet/'evidence.json').write_text(json.dumps(data))
                args+=['--packet',str(packet)]
            subprocess.run([sys.executable,str(ROOT/'py_scripts/misassembly_report.py'),*args],cwd=root,check=True,capture_output=True)
            text=(root/'README.md').read_text(encoding='utf-8')
            detail=(root/'assemblies/sample_hap1/report.md').read_text(encoding='utf-8')
            self.assertIn('sample_hap2 C01',text)
            self.assertIn('identical junctions have not been established',text)
            self.assertIn('../sample_hap2/report.md#c01',detail)
            self.assertIn('not established as the transition location',text)
            self.assertNotIn('library disagreement',text)
            self.assertEqual(table(root/'review.tsv')[0]['event_assessment'],'REVIEW')


if __name__=='__main__':unittest.main()
