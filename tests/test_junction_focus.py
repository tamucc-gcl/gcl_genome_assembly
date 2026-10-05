import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from junction_focus import read_blocks, focal_query, intersection_length, anchor_intervals, anchor_pair


class FocusTests(unittest.TestCase):
    def test_clip_and_reverse_molecule_coordinates(self):
        length,blocks=read_blocks(100,'10H20S50M5I20M',True)
        self.assertEqual(length,105)
        self.assertEqual(blocks,[(100,150,25,75),(150,170,0,20)])
        self.assertEqual(focal_query(blocks,110,120,True),[(55,65)])

    def test_deletion_not_counted_as_aligned_query(self):
        _,blocks=read_blocks(100,'20M100D20M')
        self.assertEqual(focal_query(blocks,125,200),[])
        self.assertEqual(focal_query(blocks,110,230),[(10,30)])

    def test_overlap_does_not_double_count_secondary_blocks(self):
        self.assertEqual(intersection_length([(0,10),(5,20)],[(8,12),(10,15)]),7)
        self.assertEqual(intersection_length([(0,10)],[(10,20)]),0)

    def test_partial_end_anchors_not_silently_clipped(self):
        self.assertEqual(anchor_intervals(20,30,100,25,0),{'left':None,'right':(30,55)})
        self.assertEqual(anchor_intervals(20,30,100,10,5),{'left':(5,15),'right':(35,45)})

    def test_same_target_is_not_enough_for_order(self):
        left=dict(target='s',strand='+',target_start=100,target_end=110)
        right=dict(target='s',strand='+',target_start=120,target_end=130)
        self.assertEqual(anchor_pair(left,right),('same_target_ordered',10))
        self.assertEqual(anchor_pair(right,left),('overlap_or_reordered',-30))
        self.assertEqual(anchor_pair(dict(right,strand='-'),dict(left,strand='-')),('same_target_ordered',10))
        self.assertEqual(anchor_pair(left,dict(right,strand='-'))[0],'orientation_discordant')
        self.assertEqual(anchor_pair(left,dict(right,target='t'))[0],'different_targets')

    def test_malformed_cigar_rejected(self):
        for cigar in ('*','50Mbad','0M'):
            with self.assertRaises(ValueError):read_blocks(0,cigar)


if __name__=='__main__':unittest.main()
