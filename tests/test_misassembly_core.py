import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from misassembly_core import (peer_arms,discover,invert_paf,support_grid,gap_spanning,
    summarize_peers,recommend,matched_control_sites,contact_result,read_depth,cohort_labels,independent_labels)
from misassembly_report import decisions


class MisassemblyTests(unittest.TestCase):
    def sam(self,cigar,pos=1,mq=60,nm=1,copies=2):
        d=tempfile.TemporaryDirectory();self.addCleanup(d.cleanup)
        path=Path(d.name)/'reads.sam'
        path.write_text(''.join(f'r{i}\t0\ts\t{pos}\t{mq}\t{cigar}\t*\t0\t0\t*\t*\tNM:i:{nm}\n' for i in range(copies)))
        return path

    def test_small_indels_preserve_read_continuity(self):
        for cigar in ('1000M1I1000M','1000M1D1000M','1000=1X1000='):
            self.assertEqual(support_grid(self.sam(cigar),1000,1000)['minimum'],2)

    def test_large_indels_and_skips_do_not_support_continuity(self):
        for cigar in ('1000M100D1000M','1000M100I1000M','1000M1N1000M'):
            self.assertEqual(support_grid(self.sam(cigar,nm=100),1000,1000)['minimum'],0)

    def test_verified_gap_deletion_is_counted_separately(self):
        path=self.sam('1000M100D1000M',nm=100)
        self.assertEqual(gap_spanning(path,1000,1100),2)
        self.assertEqual(gap_spanning(path,1050,1100),0)

    def test_mq_tiers_are_measured_separately(self):
        path=self.sam('4000M',mq=20)
        self.assertEqual(support_grid(path,1500,2000)['minimum'],2)
        self.assertEqual(support_grid(path,1500,2000,mapq=30)['minimum'],0)

    def test_depth_excludes_deleted_bases(self):
        rows=read_depth(self.sam('1000M1000D1000M',nm=0),0,3000)
        self.assertEqual([r['depth'] for r in rows],[2,0,2])

    def test_minor_island_does_not_create_fusion(self):
        track=[dict(lo=0,hi=2000000,chrom='chr12'),dict(lo=2000000,hi=2100000,chrom='chr4'),dict(lo=2100000,hi=4000000,chrom='chr12')]
        events,islands=peer_arms(track)
        self.assertEqual(events,[]);self.assertEqual(len(islands),1)

    def test_composite_block_labels_come_from_cohort_cigar(self):
        with tempfile.TemporaryDirectory() as tmp:
            paf=Path(tmp)/'pair.paf'
            paf.write_text('s\t400000\t0\t201000\t+\tx\t400000\t0\t200000\t200000\t201000\t60\tcg:Z:100000=1000I100000=\n')
            catalog={
                'own':dict(id='own',sample='own',taxid='1',eligible=True,scope=['s'],chromosome_labels={}),
                'peer':dict(id='peer',sample='peer',taxid='1',eligible=True,scope=['x'],chromosome_labels={'x':'chr7'})}
            result=cohort_labels(catalog,[dict(a='own',b='peer',path=str(paf))])
            self.assertEqual([(r['lo'],r['hi']) for r in result['own']['chromosome_labels']['s']],[(0,100000),(101000,201000)])
            self.assertIn('cohort',result['own']['block_label_source'])
            self.assertEqual(independent_labels(result['own'],'peer')['s'],[])

    def test_sisters_do_not_supply_independent_detection(self):
        track={'s':[dict(lo=0,hi=2000000,chrom='chr9'),dict(lo=2100000,hi=4000000,chrom='chr7')]}
        self.assertEqual(discover([(dict(id='h2',sample='x',eligible=True),track)],'x')[0],[])
        self.assertEqual(len(discover([(dict(id='p',sample='p',eligible=True),track)],'x')[0]),1)

    def test_inverse_paf_cigar_and_strand(self):
        f='a 1000 100 304 - b 2000 500 705 200 209 60 cg:Z:100M4I5D100M'.split()
        out=invert_paf(f)
        self.assertEqual(out[0],'b');self.assertEqual(out[12],'cg:Z:100M5I4D100M')
        self.assertEqual(invert_paf(out),f)

    def test_haplotype_disagreement_is_conflicting(self):
        def t(left,right):return dict(sample='p',eligible=True,qualified=True,left=dict(chrom=left),right=dict(chrom=right))
        result=summarize_peers([t('chr9','chr7'),t('chr9','chr9')],'x',('chr9','chr7'))
        self.assertEqual(result['separate'],[]);self.assertEqual(result['conflicting'],['p'])

    def test_contact_controls_match_transition_width(self):
        controls,links=matched_control_sites({'C01':dict(scaffold='s',lo=2000000,hi=2170000)},[],{'s':10000000})
        self.assertTrue(links['C01'])
        self.assertTrue(all(c['hi']-c['lo']==170000 for c in controls.values()))

    def test_no_contacts_is_not_a_loss(self):
        result=contact_result(dict(cross=0,left_within=0,right_within=0),[])
        self.assertFalse(result['loss']);self.assertIsNone(result['ratio'])

    def test_contact_loss_with_coverage_matched_controls(self):
        result=contact_result(dict(cross=0,left_within=1000,right_within=1000),[dict(id=str(i),cross=100,left_within=1000,right_within=1000) for i in range(10)])
        self.assertTrue(result['loss']);self.assertEqual(len(result['control_ids']),10)

    def test_cut_needs_observable_flanks_and_unique_option(self):
        peers=dict(separate=['a','b'],continuous=[],conflicting=[])
        option=dict(id='G1',peers=peers,spanning=0,flank_molecules=[5,5],contacts=[dict(loss=True),dict(loss=True)])
        event=dict(peers=peers,cut_options=[option],hifi=dict(state='patchy'))
        self.assertEqual(recommend(event)[0],'CUT')
        option['flank_molecules']=[0,5]
        self.assertEqual(recommend(event)[0],'SUSPECT')
        option['flank_molecules']=[5,5];event['cut_options']=[option,dict(option,id='G2')]
        self.assertEqual(recommend(event)[0],'SUSPECT')

    def test_supported_sequence_retained_despite_chromosome_change(self):
        event=dict(peers=dict(separate=['a','b'],continuous=[],conflicting=[]),hifi=dict(state='mostly_continuous'))
        self.assertEqual(recommend(event)[0],'RETAIN')

    def test_known_gap_spanning_reads_can_support_retention(self):
        peers=dict(separate=['a','b'],continuous=[],conflicting=[])
        option=dict(lo=1000,hi=1100,spanning=2,peers=peers)
        event=dict(lo=1000,hi=1100,peers=peers,hifi=dict(state='patchy'),cut_options=[option])
        self.assertEqual(recommend(event)[0],'RETAIN')
        option['lo']=2000;option['hi']=2100
        self.assertEqual(recommend(event)[0],'SUSPECT')


if __name__=='__main__':unittest.main()
