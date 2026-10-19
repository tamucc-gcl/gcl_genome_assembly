"""General junction diagnostic tests; no CTlk locations or chromosome identities."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('join_investigation',ROOT/'scripts/comparisons/join_investigation.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


class JoinInvestigationTests(unittest.TestCase):
    def test_two_haplotypes_are_not_two_independent_individuals(self):
        obs=[('B',0,100,'chr1'),('B',0,100,'chr1')]
        self.assertEqual(mod.consensus_masks(obs,'A'),[])
        self.assertEqual(mod.consensus_masks(obs+[('C',0,100,'chr1')],'A'),[dict(lo=0,hi=100,chrom='chr1')])

    def test_focal_and_conflicting_individuals_cannot_create_unique_anchors(self):
        obs=[('A',0,100,'chr1'),('B',0,100,'chr1'),('C',0,100,'chr1'),('C',25,75,'chr2')]
        self.assertEqual(mod.consensus_masks(obs,'A'),[dict(lo=0,hi=25,chrom='chr1'),dict(lo=75,hi=100,chrom='chr1')])

    def test_same_chromosome_repeat_is_distinguished_from_unique_placement(self):
        hits=[['q','10000','0','10000','+','x','10000','0','10000','10000','10000','60','cg:Z:10000='],
              ['q','10000','0','10000','+','y','10000','0','10000','10000','10000','0','cg:Z:10000=']]
        labels={'x':'chr1','y':'chr1'}
        chrom=mod.local_tracks(hits,labels,1000,2000,'chromosome')
        unique=mod.local_tracks(hits,labels,1000,2000,'locus_unique')
        self.assertEqual(chrom,[dict(lo=1000,hi=2000,chrom='chr1',status='chromosome_assigned_locus_ambiguous')])
        self.assertIsNone(unique[0]['chrom'])

    def test_reverse_cigar_clipping_keeps_target_coordinates_correct(self):
        hit=['q','10000','0','10000','-','x','10000','0','10000','10000','10000','60','cg:Z:10000=']
        labels={'x':[dict(lo=8000,hi=9000,chrom='chr2')]}
        self.assertEqual(mod.local_tracks([hit],labels,1000,2000,'locus_unique'),
                         [dict(lo=1000,hi=2000,chrom='chr2',status='assigned')])

    def test_clipping_agrees_with_production_projection(self):
        from chimera_tracks import tracks
        hits=[['q','10000','0','10000','+','x','10000','0','10000','9990','10000','60','cg:Z:10000=']]
        self.assertEqual(mod.local_tracks(hits,{'x':'chr1'},0,10000,'locus_unique'),
                         tracks(hits,{'x':'chr1'},mode='locus_unique',min_mapq=20,minimum_target=0,compact=True))

    def test_shifted_reverse_and_competing_projections_match_production(self):
        from chimera_tracks import tracks
        import random
        rng=random.Random(39)
        labels={'x':[dict(lo=0,hi=4000,chrom='chr1'),dict(lo=4000,hi=10000,chrom='chr2')],'y':'chr1'}
        for _ in range(20):
            strand=rng.choice(('+','-'));lo=rng.randrange(0,7000);hi=lo+rng.randrange(1,3000)
            hits=[['q','10000','0','10000',strand,'x','10000','0','10000','10000','10000','60','cg:Z:10000=']]
            if rng.choice((True,False)):
                hits.append(['q','10000','0','10000','+','y','10000','0','10000','9900','10000','0','cg:Z:10000='])
            for mode in ('chromosome','locus_unique'):
                full=tracks(hits,labels,mode=mode,min_mapq=20,minimum_target=0,compact=True)
                expected=[dict(b,lo=max(lo,b['lo']),hi=min(hi,b['hi'])) for b in full if b['lo']<hi and b['hi']>lo]
                self.assertEqual(mod.local_tracks(hits,labels,lo,hi,mode),expected)

    def test_breakpoint_need_not_be_last_chromosome_label(self):
        f=['r','0','s','1','60','10000=','*','0','0','A'*10000,'*','NM:i:0']
        mask=[dict(lo=0,hi=4000,chrom='chr1'),dict(lo=4500,hi=4800,chrom='chr1'),dict(lo=5000,hi=10000,chrom='chr2')]
        self.assertEqual(mod.arm_coverage(f,4000,mask,'chr1','chr2'),(4000,5000))

    def test_large_deletion_does_not_count_as_uninterrupted_bridge(self):
        f=['r','0','s','1','60','4000=1000D4000=','*','0','0','A'*8000,'*','NM:i:1000']
        self.assertFalse(mod.span_record(f,4500,1000))

    def test_unique_anchor_reach_accumulates_discontinuous_islands(self):
        mask=[dict(lo=0,hi=4000,chrom='chr1'),dict(lo=4500,hi=4800,chrom='chr1'),dict(lo=5000,hi=10000,chrom='chr2')]
        self.assertEqual(mod.minimum_anchor_span(mask,4900,'chr1','chr2'),2700)
        self.assertIsNone(mod.minimum_anchor_span(mask,4900,'chr1','chr3'))

    def test_mask_bins_rejects_partial_and_wrong_chromosome_support(self):
        mask=[dict(lo=0,hi=7000,chrom='chr1'),dict(lo=10000,hi=20000,chrom='chr2')]
        self.assertEqual(mod.mask_bins(mask,0,2,10000,'chr1'),[])
        self.assertEqual(mod.mask_bins(mask,0,2,10000,'chr2'),[1])

    def test_multiresolution_counts_preserve_within_bin_pairs(self):
        import numpy as np
        matrix=np.array([[2,3,5,7],[3,11,13,17],[5,13,19,23],[7,17,23,29]],dtype=np.uint32)
        coarse=mod.aggregate_matrix(matrix,2)
        np.testing.assert_array_equal(coarse,[[16,42],[42,71]])
        self.assertEqual(mod.window_counts(matrix,2,2),mod.window_counts(coarse,1,1))

    def test_exact_25kb_aggregation_from_5kb_counts(self):
        import numpy as np
        m=np.zeros((10,10),dtype=np.uint32);m[4,5]=m[5,4]=7;m[1,3]=m[3,1]=5
        coarse=mod.aggregate_matrix(m,5)
        np.testing.assert_array_equal(coarse,[[5,7],[7,0]])

    def test_zero_cross_is_retained_and_zero_flank_is_explicit(self):
        import numpy as np
        c=mod.window_counts(np.array([[10,0],[0,20]]),1,1)
        self.assertEqual(c['ratio'],0)
        self.assertIsNone(mod.window_counts(np.array([[0,0],[0,20]]),1,1)['ratio'])

    def test_libraries_support_one_library_and_reject_ambiguous_identity(self):
        self.assertEqual(mod.library_id('L_read',[dict(library_id='single',qname_prefix='L_')]),'single')
        with self.assertRaises(ValueError):mod.library_id('L_read',[dict(library_id='a',qname_prefix='L'),dict(library_id='b',qname_prefix='L_')])

    def test_manifest_write_preserves_timestamp_for_resume(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'input.json';mod.save(p,{'source':'unchanged'});before=p.stat().st_mtime_ns
            mod.save(p,{'source':'unchanged'});self.assertEqual(before,p.stat().st_mtime_ns)

    def test_cached_command_parsing_preserves_json_and_spaces(self):
        import shlex
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);text=json.dumps([dict(a='a',b='b',path='some path/aln.paf')])
            (p/'.command.sh').write_text('python3 misassembly_discover.py discover --alignments '+shlex.quote(text)+' --catalog '+shlex.quote('some path/catalog.json'))
            args=mod.command_args(p,'misassembly_discover.py')
            self.assertEqual(json.loads(args['alignments'])[0]['path'],'some path/aln.paf')
            self.assertEqual(args['catalog'],'some path/catalog.json')

    def test_windows_include_all_alternatives_without_known_coordinates(self):
        e=dict(id='C99',lo=100,hi=200,alternatives=[dict(id='R',lo=500,hi=600,probe_bp=550)],cut_options=[])
        sites,lo,hi=mod.windows(e,1000,padding=50)
        self.assertEqual((lo,hi),(50,650));self.assertEqual([s['point'] for s in sites],[150,550])

    @unittest.skipUnless(importlib.util.find_spec('matplotlib'),'plotting dependency not available')
    def test_dense_background_is_batched_into_two_collections(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        mask=[dict(lo=i*100,hi=i*100+50,chrom='chr1' if i%2 else 'chr2') for i in range(10000)]
        fig,ax=plt.subplots();mod.shade_masks(ax,mask,{'chr1':'blue','chr2':'orange'},.1)
        self.assertEqual(len(ax.collections),2)
        self.assertEqual(sum(len(c.get_paths()) for c in ax.collections),10000)
        plt.close(fig)


if __name__=='__main__':unittest.main()
