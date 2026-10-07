import unittest,tempfile,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_controls import literal_gaps,nearby_hypotheses
class LiteralGapTests(unittest.TestCase):
    def test_embedded_gap_crosses_lines_and_is_discovered(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.fa';p.write_text('>scaffold_1\nAAAANNNNNN\nNNNNCCCC\n>other\nNNNNNNNNNN\n')
            gaps=literal_gaps(p)
            self.assertEqual(gaps,[dict(scaffold='scaffold_1',lo=4,hi=14),dict(scaffold='other',lo=0,hi=10)])
            self.assertEqual(len(nearby_hypotheses({'c':dict(scaffold='scaffold_1',lo=0,hi=3)},gaps)),1)
    def test_short_n_runs_are_excluded(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.fa';p.write_text('>s\nAAANNNCC\n')
            self.assertEqual(literal_gaps(p),[])
