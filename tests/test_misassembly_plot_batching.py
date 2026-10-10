import sys
import unittest
from pathlib import Path
from unittest.mock import Mock

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from misassembly_assess import draw_chromosome_track


class PlotBatchingTests(unittest.TestCase):
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


if __name__=='__main__':unittest.main()
