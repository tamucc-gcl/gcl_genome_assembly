import csv
import gzip
import io
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from pairwise_paf import filter_paf,collect,parse_paf
from plot_pairwise import read_alignments,layout,overview

LINE='q\t1000\t0\t100\t+\tr\t1000\t100\t200\t90\t100\t60\n'


class PairwiseEvidenceTests(unittest.TestCase):
    def test_identity_is_matches_over_aligned_length_and_missing_dv_is_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);filter_paf(io.StringIO(LINE),root/'a.gz',root/'a.tsv','ref','qry','asm5',5,10)
            with (root/'a.tsv').open() as h:r=next(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(float(r['mean_identity']),.9);self.assertEqual(r['mean_dv'],'NA')
            self.assertEqual(r['dv_observed_alignments'],'0')

    def test_empty_alignment_is_valid_but_identity_is_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);filter_paf(io.StringIO(''),root/'a.gz',root/'a.tsv','r','q','asm5',0,1)
            with gzip.open(root/'a.gz','rt') as h:self.assertEqual(h.read(),'')
            with (root/'a.tsv').open() as h:r=next(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(r['n_align'],'0');self.assertEqual(r['mean_identity'],'NA')

    def test_divergence_uses_only_measured_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            filter_paf(io.StringIO(LINE+LINE.rstrip()+'\tdv:f:0.2\n'),root/'a.gz',root/'a.tsv','r','q','asm5',0,1)
            with (root/'a.tsv').open() as h:r=next(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(float(r['mean_dv']),.2);self.assertEqual(r['dv_observed_alignments'],'1')

    def test_zero_mapq_filter_is_honored_and_unknown_is_not_high_confidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            lines=LINE.replace('\t60\n','\t0\n')+LINE.replace('\t60\n','\t255\n')
            filter_paf(io.StringIO(lines),root/'a.gz',root/'a.tsv','r','q','asm5',0,1)
            with (root/'a.tsv').open() as h:r=next(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(r['n_align'],'2');self.assertEqual(r['mapq_observed_alignments'],'1')
            filter_paf(io.StringIO(lines),root/'a.gz',root/'a.tsv','r','q','asm5',5,1)
            with (root/'a.tsv').open() as h:r=next(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(r['n_align'],'0')

    def test_malformed_records_fail(self):
        for line in ['not a PAF',LINE.replace('\t100\t+','\t1001\t+'),
                     LINE.rstrip()+'\tdv:f:nan',LINE.rstrip()+'\tdv:f:0.1\tdv:f:0.2']:
            with self.subTest(line=line),self.assertRaises(ValueError):parse_paf(line)

    def test_empty_summary_has_schema_and_duplicate_pairs_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);collect([],root/'empty.tsv')
            self.assertTrue((root/'empty.tsv').read_text().startswith('ref_id\tqry_id'))
            for name in ['a','b']:filter_paf(io.StringIO(LINE),root/(name+'.gz'),root/(name+'.tsv'),'r','q','asm5',0,1)
            with self.assertRaisesRegex(ValueError,'Duplicate'):collect([root/'a.tsv',root/'b.tsv'],root/'out.tsv')

    def test_plot_rejects_wrong_coordinate_dictionary(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'a.gz'
            with gzip.open(path,'wt') as out:out.write(LINE)
            with self.assertRaisesRegex(ValueError,'final assembly indexes'):read_alignments(path,{'r':999},{'q':1000})
            path.write_bytes(b'bad gzip')
            with self.assertRaises(OSError):read_alignments(path,{'r':1000},{'q':1000})

    def test_axes_include_unaligned_sequences_and_both_strands_render(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rows=[parse_paf(LINE),parse_paf(LINE.replace('\t+\t','\t-\t'))]
            ref={'r':1000,'unaligned':2000};query={'q':1000}
            self.assertEqual(layout(ref,0)[2],.003)
            for kind in ['dotplot','ribbon']:
                r=overview(rows,ref,query,root/(kind+'.png'),'ref','qry',kind,0,1,5,4)
                self.assertEqual(r['alignments'],2);self.assertEqual(r['ref_axis_bp'],3000)
                self.assertGreater((root/(kind+'.png')).stat().st_size,1000)
            r=overview([],ref,query,root/'empty.png','ref','qry','ribbon',0,1,5,4)
            self.assertEqual(r['status'],'no_alignments_passing_plot_filters')


if __name__=='__main__':unittest.main()
