import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_graph_evidence import measure_graph


class NativeGraph(unittest.TestCase):
    def assay(self, directory, graph_text, hits, spanning=0):
        root=Path(directory);graph=root/'input.gfa';graph.write_text(graph_text)
        out=root/'packet';out.mkdir()
        query=root/'query.fa';query.write_text('>candidate\nAAA\n')
        intervals={'candidate':dict(start=0,lo=1000,hi=1100,hifi=dict(spanning=spanning))}
        def run(args,output):Path(output).write_text(hits)
        return measure_graph(graph,out,query,intervals,run,1),out

    def test_native_assembled_path_alone_is_not_molecule_support(self):
        with tempfile.TemporaryDirectory() as directory:
            result,out=self.assay(directory,'S\tnative\t'+'A'*3000+'\n',
                'candidate\t2100\t0\t2100\t+\tnative\t3000\t0\t2100\t2100\t2100\t60\tcg:Z:2100M\n')
            self.assertTrue(result['candidate']['native_continuity'])
            self.assertIsNone(result['candidate']['graph_contradiction'])
            self.assertTrue((out/'native_neighborhood.gfa').exists())
            self.assertFalse((out/'native_contigs.fa').exists())

    def test_unique_separate_native_segments_without_link(self):
        with tempfile.TemporaryDirectory() as directory:
            result,out=self.assay(directory,'S\tleft\t'+'A'*1000+'\nS\tright\t'+'C'*1000+'\n',
                'candidate\t2100\t0\t1000\t+\tleft\t1000\t0\t1000\t1000\t1000\t60\tcg:Z:1000M\n'
                'candidate\t2100\t1100\t2100\t+\tright\t1000\t0\t1000\t1000\t1000\t60\tcg:Z:1000M\n')
            self.assertIs(result['candidate']['graph_contradiction'],False)
            self.assertFalse(result['candidate']['native_continuity'])
            self.assertIn('S\tleft\t*\tLN:i:1000',(out/'native_neighborhood.gfa').read_text())

    def test_missing_graph_sequence_is_unavailable_not_negative(self):
        with tempfile.TemporaryDirectory() as directory:
            result,out=self.assay(directory,'S\tnative\t*\tLN:i:3000\n','')
            self.assertIsNone(result['candidate']['graph_contradiction'])
            self.assertEqual(result['candidate']['graph_status'],'native_segment_sequences_unavailable')


if __name__=='__main__':unittest.main()
