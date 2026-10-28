import csv
import gzip
import io
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from pairwise_paf import filter_paf,collect,parse_paf,union_bp

LINE='q\t1000\t0\t100\t+\tr\t1000\t100\t200\t90\t100\t60\n'


class PairwiseEvidenceTests(unittest.TestCase):
    def test_identity_is_matches_over_aligned_length_and_missing_dv_is_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);filter_paf(io.StringIO(LINE),root/'a.gz',root/'a.tsv','ref','qry','asm5',5,10)
            with (root/'a.tsv').open() as h:r=next(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(float(r['approx_match_fraction']),.9);self.assertEqual(r['mean_dv'],'NA')
            self.assertEqual(r['dv_observed_alignments'],'0')

    def test_empty_alignment_is_valid_but_identity_is_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);filter_paf(io.StringIO(''),root/'a.gz',root/'a.tsv','r','q','asm5',0,1)
            with gzip.open(root/'a.gz','rt') as h:self.assertEqual(h.read(),'')
            with (root/'a.tsv').open() as h:r=next(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(r['n_align'],'0');self.assertEqual(r['approx_match_fraction'],'NA')

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

    def test_raw_retention_and_union_coverage_do_not_double_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'r.fai').write_text('r\t1000\t0\t1000\t1001\n')
            (root/'q.fai').write_text('q\t1000\t0\t1000\t1001\n')
            lines=LINE+LINE.replace('\t100\t200\t90','\t150\t250\t90')
            filter_paf(io.StringIO(lines),root/'kept.gz',root/'qc.tsv','r','q','asm5',5,101,
                       root/'raw.gz',root/'ret.tsv',root/'r.fai',root/'q.fai')
            with gzip.open(root/'raw.gz','rt') as h:self.assertEqual(h.read(),lines)
            with gzip.open(root/'kept.gz','rt') as h:self.assertEqual(h.read(),'')
            with (root/'ret.tsv').open() as h:rows={r['filter']:r for r in csv.DictReader(h,delimiter='\t')}
            self.assertEqual(rows['unfiltered']['ref_covered_bp'],'150')
            self.assertEqual(rows['unfiltered']['query_covered_bp'],'100')
            self.assertEqual(rows['min_101']['n_align'],'0')
            self.assertEqual(union_bp({'a':[(0,10),(5,20)],'b':[(0,10)]}),30)

    def test_alignment_dictionary_mismatch_fails_before_publication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'r.fai').write_text('r\t999\t0\t999\t1000\n')
            (root/'q.fai').write_text('q\t1000\t0\t1000\t1001\n')
            with self.assertRaisesRegex(ValueError,'final assembly indexes'):
                filter_paf(io.StringIO(LINE),root/'a.gz',root/'a.tsv','r','q','asm5',5,10,
                           ref_fai=root/'r.fai',query_fai=root/'q.fai')


if __name__=='__main__':unittest.main()
