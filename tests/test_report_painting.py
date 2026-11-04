import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('painting',Path(__file__).resolve().parents[1]/'py_scripts/report_painting.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


class PaintingTests(unittest.TestCase):
    def test_overlap_consensus_and_conflict_do_not_double_count(self):
        result=mod.partition(100,[(10,40,'chr1'),(20,50,'chr1'),(30,70,'chr2')])
        self.assertEqual(result,[(0,10,'Unassigned'),(10,30,'chr1'),(30,50,'Ambiguous'),(50,70,'chr2'),(70,100,'Unassigned')])
        self.assertEqual(sum(hi-lo for lo,hi,label in result),100)

    def test_reference_composites_are_not_invented_painting_blocks(self):
        self.assertEqual(mod.chromosome('chr12_3'),'chr12')
        self.assertIsNone(mod.chromosome('chr7_1+chr12_1'))
        self.assertIsNone(mod.chromosome('unplaced_1'))

    def test_final_dictionary_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fai=root/'a.fai';paf=root/'a.paf'
            fai.write_text('chr1_1\t100\t0\t100\t101\n')
            paf.write_text('chr1_1\t101\t0\t80\t+\tchr1_1\t100\t0\t80\t80\t80\t60\n')
            with self.assertRaisesRegex(ValueError,'final sequence dictionaries'):
                mod.paint(fai,fai,paf,root/'out.tsv',min_length=1)

    def test_mapq_255_is_unavailable_and_self_is_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fai=root/'a.fai';paf=root/'a.paf';out=root/'out.tsv'
            fai.write_text('chr1_1\t100\t0\t100\t101\n')
            paf.write_text('chr1_1\t100\t0\t80\t+\tchr1_1\t100\t0\t80\t80\t80\t255\n')
            mod.paint(fai,fai,paf,out,min_length=1)
            self.assertIn('Unassigned',out.read_text())
            result=mod.paint(fai,fai,paf,out,self_reference=True)
            self.assertTrue(result['self_reference'])
            self.assertIn('reference_identity',out.read_text())


if __name__=='__main__':unittest.main()
