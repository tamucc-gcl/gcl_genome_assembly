import unittest,sys,tempfile,json,csv,subprocess
from pathlib import Path
from collections import defaultdict,Counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_blocks import chromosome_blocks,summarize_blocks,farther_regions,farther_contacts
from chimera_adjudicate import decide

class Blocks(unittest.TestCase):
    def measurement(self):
        return dict(verified_gap=True,review_only=True,gap_start=10000000,gap_end=10000100,
            cut_bp=10000100,assessment_scaffold_length=20000100,assessment_sha256='a'*64,
            hifi_spanning_molecules=0,hifi_informative=False,alternative_placements_checked=True,
            graph_contradiction=False,haplotype_block_conflict=False,repeat_obscured_localization=True,
            chromosome_blocks=dict(localized=True,independent_individuals=['a','b','c']),
            farther_contact_evidence=dict(pass_all=True))
    def test_obscured_candidate_can_qualify_without_immediate_flank_measurements(self):
        self.assertEqual(decide(self.measurement())[:2],('UNJOIN_UNSUPPORTED',True))
    def test_intact_and_ambiguous_and_conflicting_cases_do_not_cut(self):
        for change in [dict(hifi_spanning_molecules=2),dict(graph_contradiction=True),
            dict(repeat_obscured_localization=False),dict(haplotype_block_conflict=True),
            dict(farther_contact_evidence=dict(pass_all=False))]:
            self.assertFalse(decide(dict(self.measurement(),**change))[1])
    def test_constructed_two_chromosome_join_localizes_only_one_gap(self):
        interval=dict(scaffold='s',lo=1000000,hi=1000100,start=150000,end=1850100)
        # Deliberately constructed cross-chromosome query, with independently placeable sides.
        hits=[['q','1700100','0','850000','+','one','20000000','0','850000','850000','850000','60','cg:Z:850000M'],
              ['q','1700100','850100','1700100','+','two','20000000','0','850000','850000','850000','60','cg:Z:850000M']]
        peer=dict(id='p',sample='a',auto_evidence=True,chromosome_labels={'one':'chr1','two':'chr2'})
        gap=dict(scaffold='s',lo=1000000,hi=1000100)
        row=chromosome_blocks(hits,interval,peer,[gap])
        self.assertTrue(row['qualified']);self.assertTrue(row['unique_gap_localization'])
        rows=[dict(row,sample=s) for s in ['a','b','c']]
        self.assertTrue(summarize_blocks(rows,'target')['localized'])
        self.assertFalse(chromosome_blocks(hits,interval,peer,[gap,dict(scaffold='s',lo=950000,hi=950100)])['unique_gap_localization'])
        continuous=dict(row,sample='d',chromosome_pair=['chr1','chr1'])
        self.assertFalse(summarize_blocks(rows+[continuous],'target')['localized'])
    def test_farther_controls_require_both_libraries_and_multiple_offsets(self):
        h=dict(informative=True)
        candidate=dict(key='x',scaffold='s',lo=1000000,hi=1000100,farther_hifi={str(o):h for o in [100000,250000,500000]})
        controls=[]
        for role in ['continuous_control','gap_control']:
            for i in range(5):controls.append(dict(candidate,key=role+str(i),role=role,hifi=dict(spanning=2),peer_continuous_individuals=3))
        intervals={c['key']:c for c in [candidate]+controls};regions=farther_regions(intervals,{'s':4000000})
        counts=defaultdict(Counter)
        for lib in ['A','B']:
            for key in regions:counts[lib,key].update(left_within=1000,right_within=1000,cross=0 if key.startswith('x@') else 100)
        self.assertTrue(farther_contacts(candidate,controls,counts,{'A','B'},regions)['pass_all'])
        counts['B','x@250000']['cross']=100
        self.assertFalse(farther_contacts(candidate,controls,counts,{'A','B'},regions)['pass_all'])
    def test_hypothesis_action_reaches_source_bound_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);context=root/'context';context.mkdir();calls=root/'calls.tsv'
            calls.write_text('scaffold\tassembly_sha256\n')
            m=dict(self.measurement(),assembly='sample',scaffold='s',coordinate_stage='pre_finishing',packet_interval_id='hyp')
            (context/'decision_measurements.json').write_text(json.dumps({'hypothesis':m}))
            script=Path(__file__).resolve().parents[1]/'py_scripts/chimera_adjudicate.py'
            subprocess.run([sys.executable,str(script),'--calls',str(calls),'--assembly','sample','--context',str(context),'--out',str(root/'out')],check=True)
            with (root/'out/actions.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(len(rows),1);self.assertEqual(rows[0]['cut_bp'],'10000100');self.assertEqual(rows[0]['policy_version'],'gap-v2')
