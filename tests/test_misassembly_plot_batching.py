import sys
import unittest
from pathlib import Path
from unittest.mock import Mock,patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from misassembly_assess import contact_profile,assessment_gaps
from misassembly_plot import draw_chromosome_track,shade_chromosomes,chromosome_colors,plot_event,plot_scaffold,detail_window
from misassembly_core import consensus_regions,discover,recommend
import tempfile


class PlotBatchingTests(unittest.TestCase):
    def test_closeup_prioritizes_user_coordinate_then_recommended_gap(self):
        e=dict(lo=1000000,hi=9000000,recommended_cut=dict(lo=3999950,hi=4000050,cut_bp=4000050))
        self.assertEqual(detail_window(e),(3900050,4100050))
        self.assertEqual(detail_window(e,[dict(cut_bp=6000000)]),(5900000,6100000))

    def test_batching_preserves_every_interval_and_unassigned_gap(self):
        axis=Mock();colors={'chr1':'red','chr2':'blue'}
        blocks=[dict(lo=0,hi=100,chrom='chr1'),dict(lo=100,hi=200,chrom=None),
                dict(lo=200,hi=300,chrom='chr2'),dict(lo=300,hi=400,chrom='chr1')]
        self.assertEqual(draw_chromosome_track(axis,blocks,2,colors),2)
        calls=axis.broken_barh.call_args_list
        self.assertEqual(calls[0].args,([(0.,.0001),(.0003,.0001)],(1.65,.7)))
        self.assertEqual(calls[1].args,([(.0002,.0001)],(1.65,.7)))
        self.assertEqual(calls[0].kwargs['facecolors'],'red')
        self.assertEqual(calls[1].kwargs['facecolors'],'blue')

    def test_collection_count_depends_on_chromosomes_not_interval_count(self):
        axis=Mock();blocks=[dict(lo=3*i,hi=3*i+2,chrom='chr1') for i in range(100000)]
        self.assertEqual(draw_chromosome_track(axis,blocks,0,{'chr1':'red'}),1)
        axis.broken_barh.assert_called_once()
        self.assertEqual(len(axis.broken_barh.call_args.args[0]),100000)

    def test_consensus_excludes_sister_and_conflicting_haplotypes(self):
        def t(peer,sample,chrom):return dict(peer=peer,sample=sample,eligible=True,blocks=[dict(lo=0,hi=10000,chrom=chrom)])
        tracks=[t('a1','a','chr1'),t('a2','a','chr2'),t('b','b','chr1'),t('own2','own','chr1')]
        self.assertEqual(consensus_regions(tracks,'own'),[])
        tracks += [t('c','c','chr1')]
        self.assertEqual(consensus_regions(tracks,'own'),[dict(lo=0,hi=10000,chrom='chr1',aligned_bp=10000)])

    def test_background_leaves_ambiguous_interval_unshaded(self):
        axis=Mock()
        shade_chromosomes([axis],[dict(lo=0,hi=100,chrom='chr1'),dict(lo=200,hi=300,chrom='chr2')],{'chr1':'red','chr2':'blue'})
        self.assertEqual([c.args for c in axis.axvspan.call_args_list],[(0,.0001),(.0002,.0003)])

    def test_colors_stable_across_chromosome_subsets_and_unique(self):
        colors=chromosome_colors(['chr'+str(i) for i in range(1,25)])
        self.assertEqual(colors['chr14'],chromosome_colors(['chr14'])['chr14'])
        self.assertEqual(len(set(colors.values())),24)

    def test_contact_profile_missing_coverage_is_unavailable_not_zero(self):
        import numpy as np
        matrix=np.zeros((6,6),dtype=np.uint32)
        self.assertIsNone(contact_profile(matrix)[0]['ratio'])
        matrix[:3,:3]=10;matrix[3:,3:]=10
        row=contact_profile(matrix)[0]
        self.assertEqual(row['ratio'],0)
        self.assertEqual(row['left_within'],60)
        matrix[:3,3:]=2;matrix[3:,:3]=2
        self.assertEqual(contact_profile(matrix)[0]['ratio'],.3)

    def test_literal_gaps_include_inherited_runs_and_wrapped_current_gap(self):
        with tempfile.TemporaryDirectory() as tmp:
            fa=Path(tmp)/'a.fa';fa.write_text('>s\nAAAAANNNNN\nNNNNNAAAAANNNNNNNNNNAAAAA\n')
            gaps=assessment_gaps(fa,[dict(scaffold='s',lo=5,hi=15)])
            self.assertEqual([(g['lo'],g['hi'],g['origin']) for g in gaps],[(5,15,'current_agp'),(20,30,'inherited_sequence_gap')])
            with self.assertRaises(ValueError):assessment_gaps(fa,[dict(scaffold='s',lo=4,hi=15)])

    def test_regional_discovery_ignores_repeated_minor_label_islands(self):
        blocks=[]
        for i in range(80):
            lo=i*1000;chrom='chr1' if i<40 else 'chr2'
            blocks += [dict(lo=lo,hi=lo+990,chrom=chrom),dict(lo=lo+990,hi=lo+1000,chrom='chr3')]
        peers=[(dict(id=p,sample=p,eligible=True),{'s':blocks}) for p in ('a','b')]
        events,_=discover(peers,'own',minimum=10000)
        self.assertEqual(len(events),1)
        self.assertEqual((events[0]['left'],events[0]['right']),('chr1','chr2'))
        self.assertLessEqual(events[0]['hi']-events[0]['lo'],1000)

    def test_single_peer_transition_opposed_by_cohort_is_audited(self):
        transition=[dict(lo=0,hi=40000,chrom='chr1'),dict(lo=40000,hi=80000,chrom='chr2')]
        continuous=[dict(lo=0,hi=80000,chrom='chr1')]
        peers=[(dict(id=p,sample=p,eligible=True),{'s':transition if p=='a' else continuous}) for p in ('a','b','c')]
        events,audit=discover(peers,'own',minimum=10000)
        self.assertEqual(events,[])
        self.assertTrue(any(r.get('kind')=='rejected_peer_transition' for r in audit))

    def test_main_figure_has_all_evidence_and_contact_map_stays_separate(self):
        import matplotlib.pyplot as plt
        import numpy as np
        event=dict(id='C01',scaffold='s',lo=400000,hi=401000,left='chr1',right='chr2',peer_tracks=[],cut_options=[],
                   contact_profiles={'A':[dict(position=400000,ratio=.2)],'B':[dict(position=400000,ratio=0)]},
                   telomeres=[dict(start=400000,count=5)])
        figures=[]
        def capture(fig,*args,**kwargs):figures.append(fig)
        with tempfile.TemporaryDirectory() as tmp,patch('matplotlib.figure.Figure.savefig',capture):
            out=Path(tmp);(out/'plots').mkdir()
            plot_event(event,out)
            self.assertEqual(len(figures[0].axes),5)
            self.assertEqual([line.get_label() for line in figures[0].axes[3].lines],['A','B'])
            self.assertEqual(figures[0].axes[4].lines[0].get_ydata().tolist(),[5])
            plot_scaffold('s',np.eye(8,dtype=np.uint32),[event],out)
            self.assertEqual(len(figures[1].axes),2) # Matrix axis and its colorbar only.
        self.assertFalse(plt.get_fignums())

    def test_region_and_selected_cut_markers_cross_every_panel(self):
        event=dict(id='C01',scaffold='s',lo=400000,hi=401000,left='chr1',right='chr2',peer_tracks=[],cut_options=[],
            alternatives=[dict(id='C01R01',lo=399000,hi=401000,probe_bp=400000,site_verdict='REVIEW')])
        figures=[]
        with tempfile.TemporaryDirectory() as tmp,patch('matplotlib.figure.Figure.savefig',lambda fig,*a,**kw:figures.append(fig)):
            out=Path(tmp);(out/'plots').mkdir()
            plot_event(event,out,detail=True,selected_cuts=[dict(cut_bp=402000)])
        for axis in figures[0].axes:
            self.assertTrue(any(line.get_color()=='#762a83' and line.get_linestyle()=='--' for line in axis.lines))
            self.assertTrue(any(line.get_color()=='#007A3D' for line in axis.lines))
        self.assertIn('Hi-C contact support',figures[0].axes[3].get_ylabel())
        self.assertEqual(figures[0].axes[4].get_xlabel(),'Original scaffold coordinate (Mb)')

    def test_ambiguous_spike_does_not_set_primary_read_scale(self):
        event=dict(assembly='assembly_hap1',id='C01',scaffold='s',lo=400000,hi=401000,left='chr1',right='chr2',peer_tracks=[],cut_options=[],
                   hifi=dict(probes=[dict(position=400000,molecules=10),dict(position=401000,molecules=8)]),
                   hifi_all_primary=dict(probes=[dict(position=400000,molecules=400),dict(position=401000,molecules=8)]),
                   depth=[dict(start=400000,end=401000,depth=10)],depth_all_primary=[dict(start=400000,end=401000,depth=400)])
        figures=[]
        with tempfile.TemporaryDirectory() as tmp,patch('matplotlib.figure.Figure.savefig',lambda fig,*a,**kw:figures.append(fig)):
            out=Path(tmp);(out/'plots').mkdir();plot_event(event,out)
        self.assertLess(figures[0].axes[1].get_ylim()[1],20)
        self.assertGreater(figures[0].axes[5].get_ylim()[1],400)
        self.assertLess(figures[0].axes[2].get_ylim()[1],20)
        self.assertIn('assembly_hap1',figures[0]._suptitle.get_text())

    def test_supported_regional_read_scan_does_not_hide_unsupported_nearby_gap(self):
        peers=dict(separate=['a','b'],continuous=[],conflicting=[])
        gap=dict(id='G1',lo=1010,hi=1020,spanning=0,flank_molecules=[5,5],peers=peers,flank_support=[dict(observable=True),dict(observable=True)],association=dict(status='confirmed'),contacts=[])
        event=dict(lo=1000,hi=1005,peers=peers,cut_options=[gap],hifi=dict(state='continuous'),hifi_strict=dict(state='continuous'))
        self.assertEqual(recommend(event)[0],'SUSPECT')
        gap['flank_support'][0]['observable']=False
        self.assertEqual(recommend(event)[0],'RETAIN')


if __name__=='__main__':unittest.main()
