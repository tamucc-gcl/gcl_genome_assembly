"""Known-answer interval/filter checks for the removable Dipcall investigation."""
import gzip
import importlib.util
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('psmc_mask_investigation',Path(__file__).resolve().parents[1]/'scripts/comparisons/psmc_mask_investigation.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
try:
    import numpy
    import pysam
    READERS=True
except ImportError:READERS=False


def paf(name,lo,hi,length=None,mapq=60,tp='P',scoring='s1:i:100\ts2:i:0'):
    return f'{name}\t300000\t0\t{hi-lo}\t+\tchr1\t300000\t{lo}\t{hi}\t{hi-lo}\t{length or hi-lo}\t{mapq}\ttp:A:{tp}\t{scoring}\tcs:Z::{hi-lo}'


class Investigation(unittest.TestCase):
    def test_filter_boundaries_use_paf_block_length(self):
        for length,quality,expected in [(9999,60,2),(10000,4,3),(10000,5,4),(49999,60,4),(50000,5,5)]:
            self.assertEqual(mod.paf_level(paf('a',0,100000,length,quality).split('\t')),expected)
        self.assertEqual(mod.paf_level(paf('a',0,100000,tp='S').split('\t')),0)
        self.assertEqual(mod.paf_level(paf('a',0,100000,tp='i').split('\t')),0)
        self.assertEqual(mod.paf_level(paf('a',0,100000,scoring='s1:i:100').split('\t')),1)
        self.assertEqual(mod.paf_level(paf('a',0,100000,tp='I').split('\t')),5)

    def test_overlap_distinguishes_one_query_from_multiple_queries(self):
        for query,expected in [('a','overlapping_segments_one_query'),('b','overlapping_multiple_queries')]:
            events=[(0,1,-1,''),(100,-1,-1,''),(0,1,5,'a'),(70,-1,5,'a'),
                    (30,1,4,query),(100,-1,4,query)]
            pieces=list(mod.sweep(events))
            self.assertEqual([(lo,hi,state) for lo,hi,state,_ in pieces],
                [(0,30,'unique_long_alignment'),(30,70,expected),(70,100,'unique_10_to_50kb')])

    def test_missing_coverage_and_scope_gaps_are_not_ambiguous_coverage(self):
        events=[(0,1,-1,''),(10,-1,-1,''),(20,1,-1,''),(40,-1,-1,''),
                (25,1,5,'a'),(40,-1,5,'a'),(30,1,-2,''),(40,-1,-2,'')]
        pieces=list(mod.sweep(events))
        self.assertEqual(sum(hi-lo for lo,hi,_,_ in pieces),30)
        self.assertEqual(sum(hi-lo for lo,hi,state,_ in pieces if state=='no_PAF_coverage'),15)
        self.assertEqual(sum(hi-lo for lo,hi,_,official in pieces if official),10)

    def test_short_or_low_quality_overlap_does_not_remove_eligible_long_alignment(self):
        events=[(0,1,-1,''),(100,-1,-1,''),(0,1,5,'a'),(100,-1,5,'a'),
                (20,1,2,'b'),(50,-1,2,'b'),(30,1,3,'c'),(80,-1,3,'c')]
        self.assertTrue(all(state=='unique_long_alignment' for _,_,state,_ in mod.sweep(events)))

    def test_interval_intersection_has_no_neighbour_filling(self):
        self.assertEqual(mod.intersect([(0,10),(20,30)],[(5,25)]),[(5,10),(20,25)])
        self.assertEqual(mod.merge([(20,25),(0,10),(10,15)]),[(0,15),(20,25)])

    @unittest.skipUnless(READERS and shutil.which('paftools.js'),'Requires pinned native readers and paftools')
    def test_predicted_coverage_matches_native_paftools_mask(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fa=root/'ref.fa';fa.write_text('>chr1\n'+'A'*300000+'\n');pysam.faidx(str(fa))
            lines=[paf('a',0,100000),paf('b',50000,150000),paf('c',150000,190000),
                paf('d',190000,200000),paf('e',200000,210000,mapq=4),paf('f',210000,270000,tp='S'),
                paf('g',270000,290000,scoring='s1:i:100')]
            plain=root/'alignment.paf';plain.write_text('\n'.join(lines)+'\n')
            run=subprocess.run(['paftools.js','call',str(plain)],capture_output=True,text=True,check=True)
            mask=root/'caller.bed';mask.write_text(''.join('\t'.join(line.split('\t')[1:])+'\n' for line in run.stdout.splitlines() if line.startswith('R\t')))
            compressed=root/'alignment.paf.gz'
            with gzip.open(compressed,'wt') as f:f.write(plain.read_text())
            job=dict(filters=dict(min_contig_bp=1),sample='test',haplotype='hap2',reference=dict(sha256='a'*64),expected_scope_bp=300000)
            with pysam.FastaFile(str(fa)) as reader:r=mod.analyze(compressed,mask,reader,job,root)
            self.assertEqual(r['official_only_bp']+r['predicted_only_bp'],0)
            self.assertEqual(r['official_confident_bp'],100000)
            self.assertEqual(r['classification_bp']['overlapping_multiple_queries'],50000)
            self.assertEqual(sum(r['classification_bp'].values()),300000)


if __name__=='__main__':unittest.main()
