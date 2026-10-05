"""Decision safeguards and complete-inventory regression tests (no external tools)."""
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'py_scripts'))
sys.path.insert(0,str(ROOT/'scripts/comparisons'))
from junction_review import n_runs, unique_anchor, compare_anchors, recommend, apply_reviews
from review_all_junctions import inventory, reuse_packets, specifications


def hit(target='s',start=0,strand='+',mapq=60):
    return dict(target=target,target_start=start,target_end=start+1000,
                strand=strand,mapq=mapq,query_length=1000,cigar='1000=')


class ReviewTests(unittest.TestCase):
    def test_ambiguous_and_short_hits_are_not_unique(self):
        self.assertIsNone(unique_anchor([hit(),hit('repeat')]))
        self.assertIsNone(unique_anchor([hit(mapq=255)]))
        self.assertIsNone(unique_anchor([dict(hit(),cigar='200=800I')]))
        self.assertIsNotNone(unique_anchor([hit()]))

    def test_orientation_and_target_comparison(self):
        self.assertEqual(compare_anchors([hit(start=2000,strand='-')],[hit(strand='-')])['status'],'same_target_ordered')
        self.assertEqual(compare_anchors([hit()],[hit('other')])['status'],'different_targets')
        self.assertEqual(compare_anchors([hit()],[hit(start=2000,strand='-')])['status'],'orientation_discordant')

    def test_sister_conflict_never_automatically_breaks(self):
        sister=[dict(status='different_targets',anchor_context='clean')]
        self.assertEqual(recommend({},sister)[0],'PRIORITIZE_MISJOIN_REVIEW')
        self.assertEqual(recommend({'hifi_chain_molecules':2},sister)[0],'REVIEW_CONFLICT')
        sister[0]['anchor_context']='confounded'
        self.assertEqual(recommend({},sister)[0],'NOT_RESOLVED')

    def test_zero_support_and_missing_not_cut_evidence(self):
        for value in ('.',0):
            self.assertEqual(recommend({'hifi_chain_molecules':value},[])[0],'NOT_RESOLVED')

    def test_review_is_bound_to_sequence_and_safe_boundary(self):
        row=dict(id='g',assembly='a',scaffold='s',start=10,end=20,
                 assessment_sha256='abc',kind='verified_agp_gap')
        review=dict(row,decision='UNJOIN_UNSUPPORTED',reviewer='person',reason='reason',
                    evidence_for='reviewed',evidence_against='none observed')
        apply_reviews([row],[review])
        self.assertEqual(row['decision'],'UNJOIN_UNSUPPORTED')
        for change in ({'assessment_sha256':'wrong'},{'start':11},{'reason':''}):
            with self.assertRaises(ValueError): apply_reviews([row],[dict(review,**change)])
        with self.assertRaises(ValueError):
            apply_reviews([dict(row,kind='flagged_interval')],[review])

    def test_inventory_includes_unflagged_scaffold_and_unknown_N_run(self):
        class FakeFasta:
            lengths={'s1':12,'s2':10,'s3':4}
            sequences={'s1':'AAAANNAAAAAA','s2':'AANNNAAAAA','s3':'AAAA'}
            def __init__(self,*args): pass
            def fetch(self,name,lo,hi): return self.sequences[name][lo:hi]
        self.assertEqual(n_runs('aaNNnAaN'),[(2,5),(7,8)])
        with tempfile.TemporaryDirectory() as tmp:
            a=SimpleNamespace(out=Path(tmp),registry=None,samtools='unused')
            specs=[dict(assembly='a',sample='fish',haplotype='1',sister='',
                        stages=[dict(stage='assessment',fasta='unused',agp='')])]
            gap=dict(id='ignored',assembly='a',scaffold='s1',start=4,end=6,
                     origin_stage='round1',left_component='c1',right_component='c2')
            with patch('review_all_junctions.sha',return_value='hash'),patch('review_all_junctions.Fasta',FakeFasta),patch('review_all_junctions.find_gaps',return_value=[gap]):
                jobs,rows=inventory(a,specs)
            self.assertEqual(len(rows),2)
            self.assertEqual({r['kind'] for r in rows},{'verified_agp_gap','unresolved_N_run'})
            self.assertIn('s3',(Path(tmp)/'scaffold_inventory.tsv').read_text())
            self.assertTrue(all(r['decision']=='UNREVIEWED' for r in rows))

    def test_sister_cannot_be_another_individual(self):
        import json
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'manifest.json'
            path.write_text(json.dumps([dict(assembly='a',sample='x',sister='b',stages=[dict(stage='assessment')]),dict(assembly='b',sample='y',sister='a',stages=[dict(stage='assessment')])]))
            with self.assertRaises(ValueError): specifications(SimpleNamespace(manifest=path))

    def test_stale_read_packet_rejected(self):
        import json
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp); (path/'a').mkdir()
            (path/'status.json').write_text(json.dumps(dict(status='SUCCESS',assemblies=[dict(assembly='a',folder='a',assessment_sha256='old')])))
            (path/'a/hifi_support.tsv').write_text('id\tscaffold\tstart\tend\tmolecules_passing_chain_screen\ng\ts\t1\t2\t3\n')
            row=dict(id='g',assembly='a',assessment_sha256='new',scaffold='s',start=1,end=2)
            with self.assertRaises(ValueError): reuse_packets(SimpleNamespace(evidence_packet=[path]),[row])


if __name__=='__main__': unittest.main()
