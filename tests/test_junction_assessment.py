"""Synthetic coordinate/evidence checks. No real genomes or installed aligners required."""
import json
import gzip
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'py_scripts'))
sys.path.insert(0,str(ROOT/'scripts/comparisons'))
from junction_assessment import (validate_registry,stable_id,project_gap,choose_gaps,Lift,
                                 library_for,Contacts,decision_rows,validate_decisions)
from assess_junction_batch import find_gaps,bam_frame,pairs_dictionary
from chimera_origin import agp_rows


def seed(id='J1',start=100,end=200):
    return dict(id=id,sample='fish',assembly='fish_hap1',scaffold='s',start=str(start),end=str(end),
                kind='assignment_transition',assessment_sha256='a'*64)


def gap(id='g1',start=140,end=150):
    return dict(id=id,assembly='fish_hap1',scaffold='s',start=start,end=end)


class JunctionTests(unittest.TestCase):
    def test_pairs_dictionary_requires_unique_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'pairs.gz'
            with gzip.open(path,'wt') as out:
                out.write('#chromsize: s 100\nread\ts\t1\ts\t2\t+\t-\tUU\n')
            self.assertEqual(pairs_dictionary(path),{'s':100})
            for content in ('#chromsize: s 100\n#chromsize: s 100\n','read\ts\t1\ts\t2\t+\t-\tUU\n'):
                with gzip.open(path,'wt') as out:
                    out.write(content)
                with self.assertRaises(ValueError):
                    pairs_dictionary(path)

    def test_registry_duplicate_and_wrong_hash(self):
        self.assertEqual(len(validate_registry([seed()])),1)
        for rows in ([seed(),seed()],[dict(seed(),assessment_sha256='bad')],[dict(seed(),assembly='../outside')]):
            with self.assertRaises(ValueError):
                validate_registry(rows)

    def test_reverse_gap_projection(self):
        self.assertEqual(project_gap(dict(start=20,end=30),100,True),(70,80))

    def test_identity_includes_assembly_and_interval(self):
        self.assertNotEqual(stable_id('a','s',1,2),stable_id('b','s',1,2))
        self.assertEqual(stable_id('a','s',1,2),stable_id('a','s',1,2))

    def test_no_nearest_gap_snapping(self):
        selected,links=choose_gaps([seed()],[gap(start=201,end=202)])
        self.assertEqual(selected[0]['role'],'nearby_control')
        self.assertTrue(any(l['relationship']=='no_verified_gap_in_interval' for l in links))

    def test_candidate_cannot_be_control_for_another_seed(self):
        selected,links=choose_gaps([seed(),seed('J2',300,400)],[gap(),gap('g2',340,350)])
        self.assertTrue(all(g['role']=='candidate_gap' for g in selected))
        self.assertFalse(any(l['relationship'].startswith('nearby_control') for l in links))

    def test_endpoint_touch_is_not_overlap(self):
        selected,links=choose_gaps([seed()],[gap(start=200,end=210)],controls=0)
        self.assertEqual(selected,[])
        self.assertEqual(links[0]['relationship'],'no_verified_gap_in_interval')

    def test_pairs_reverse_component_and_duplicate_placement(self):
        r=dict(kind='W',component='c',component_start=10,component_end=30,object='s',start=100,end=120,orientation='-')
        self.assertEqual(Lift([r]).locate('c',10),('s',119))
        self.assertEqual(Lift([r]).locate('c',29),('s',100))
        self.assertIsNone(Lift([r]).locate('c',30))
        self.assertIsNone(Lift([r,dict(r,object='t')]).locate('c',15))

    def test_library_groups_readsets_without_merging_libraries(self):
        manifest=[dict(qname_prefix='r1_',library_id='A'),dict(qname_prefix='r2_',library_id='A'),dict(qname_prefix='r3_',library_id='B')]
        self.assertEqual(library_for('r2_read',manifest),'A')
        self.assertEqual(library_for('r3_read',manifest),'B')
        self.assertEqual(library_for('oldname',manifest),'UNASSIGNED')
        with self.assertRaises(ValueError):
            library_for('r1_read',manifest+[dict(qname_prefix='r',library_id='C')])

    def test_contacts_pair_order_zero_and_library_separation(self):
        a=Contacts([gap(start=100,end=110)],{'s':300,'t':300},50)
        a.add('A',('s',60),('s',120))
        a.add('A',('s',120),('s',60))
        a.add('B',('s',60),('t',10))
        self.assertEqual(a.cross['A','g1'],2)
        self.assertEqual(a.cross['B','g1'],0)
        self.assertEqual(a.partners['B','g1','left','t',0],1)
        a.add('A',('s',60),('s',70))
        self.assertEqual(a.margin['A','g1','left'],2)

    def test_contacts_half_open_and_clipped_flanks(self):
        a=Contacts([gap(start=10,end=20)],{'s':35},50)
        self.assertEqual(a.windows['g1','left'],('s',0,10))
        self.assertFalse(a.ends('s',10))
        self.assertEqual(a.ends('s',20),{('g1','right')})

    def test_decisions_do_not_authorize_cuts(self):
        rows=decision_rows([seed()],[dict(seed='J1',gap='g1',relationship='intersects_seed_interval')])
        self.assertEqual(rows[0]['cut_authorized'],'no')
        self.assertTrue(validate_decisions(rows,[seed()],[gap()]))
        rows[0].update(action='BREAK',chosen_gap='g1')
        with self.assertRaises(ValueError):
            validate_decisions(rows,[seed()],[gap()])
        rows[0].update(reviewer='reviewer',rationale='reviewed evidence',evidence_for='packet',evidence_against='none observed')
        self.assertTrue(validate_decisions(rows,[seed()],[gap()]))
        with self.assertRaises(ValueError):
            validate_decisions(rows,[seed()],[gap(start=500,end=510)])

    def test_wrong_bam_provenance_rejected_before_external_commands(self):
        with patch('assess_junction_batch.run') as run:
            with self.assertRaises(ValueError):
                bam_frame(dict(assembly='wrong'),'fish','a'*64,Path('bam'),None,'samtools')
            run.assert_not_called()

    def test_literal_stage_projection_and_changed_stage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            fa=root/'stage.fa';fa.write_text('placeholder')
            agp=root/'stage.agp'
            agp.write_text('s\t1\t4\t1\tW\tc1\t1\t4\t+\ns\t5\t6\t2\tN\t2\tscaffold\tyes\tproximity_ligation\ns\t7\t10\t3\tW\tc2\t1\t4\t+\n')
            class FakeFasta:
                lengths={'s':10}
                def __init__(self,seq):self.seq=seq
                def fetch(self,name,lo,hi):return self.seq[lo:hi]
            specs=[dict(stage='round1',fasta=str(fa),agp=str(agp)),dict(stage='assessment',fasta=str(fa),agp=str(agp))]
            with patch('assess_junction_batch.Fasta',return_value=FakeFasta('AAAANNGGGG')):
                gaps=find_gaps('fish_hap1',FakeFasta('AAAANNGGGG'),specs,root,'samtools',{'s'})
            self.assertEqual(len(gaps),1)
            self.assertEqual((gaps[0]['start'],gaps[0]['end'],gaps[0]['origin_stage']),(4,6,'round1'))
            with patch('assess_junction_batch.Fasta',return_value=FakeFasta('TTTTNNGGGG')):
                gaps=find_gaps('fish_hap1',FakeFasta('AAAANNGGGG'),specs,root,'samtools',{'s'})
            self.assertEqual(gaps[0]['origin_stage'],'assessment')
            with patch('assess_junction_batch.Fasta',return_value=FakeFasta('AAAACCGGGG')):
                with self.assertRaises(ValueError):
                    find_gaps('fish_hap1',FakeFasta('AAAACCGGGG'),specs,root,'samtools',{'s'})


if __name__=='__main__':
    unittest.main()
