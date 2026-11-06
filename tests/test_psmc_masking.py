"""Known-genotype/unknown-site assertions with actual pysam VCF/FASTA readers."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from psmc_input import prepare,finalize
try:
    import pysam
    import numpy
    AVAILABLE=True
except ImportError:
    AVAILABLE=False


@unittest.skipUnless(AVAILABLE,'Native input checks require pinned pysam/numpy environment')
class Masking(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.fa=self.root/'reference.fa'
        self.fa.write_text('>long\n'+'A'*200+'N'*20+'A'*180+'\n>short\n'+'A'*10+'\n')
        pysam.faidx(str(self.fa))
        self.args=argparse.Namespace(reference=str(self.fa),reference_sha256=hashlib.sha256(self.fa.read_bytes()).hexdigest(),
            mode='assemblies',min_contig_bp=100,min_depth=4,max_depth=20,min_gq=20,indel_padding=5,
            min_allele_balance=.2,include=None,exclude=None,depth=None,confident=None,out=str(self.root/'out'))
    def tearDown(self):self.temp.cleanup()
    def vcf(self,rows,formats):
        header=('##fileformat=VCFv4.2\n##contig=<ID=long,length=400>\n'
                '##FILTER=<ID=LowQual,Description="low">\n'+formats+
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n')
        p=self.root/'calls.vcf';p.write_text(header+''.join('long\t'+r+'\n' for r in rows))
        self.args.vcf=pysam.tabix_index(str(p),preset='vcf',force=True)
    def consensus(self):
        with gzip.open(self.root/'out/diploid.fa.gz','rt') as f:
            return ''.join(l.strip() for l in f if not l.startswith('>'))
    def test_dipcall_scope_gaps_failed_genotypes_and_complex_flanks(self):
        conf=self.root/'confident.bed';conf.write_text('long\t10\t390\n');self.args.confident=str(conf)
        exclude=self.root/'exclude.bed';exclude.write_text('long\t100\t110\n');self.args.exclude=str(exclude)
        self.vcf(['21\t.\tA\tG\t60\tPASS\t.\tGT\t0/1',
                  '51\t.\tAA\tA\t60\tPASS\t.\tGT\t0/1',
                  '80\t.\tA\tG\t60\tPASS\t.\tGT\t0/.',
                  '150\t.\tA\tG\t10\tLowQual\t.\tGT\t0/1',
                  '300\t.\tA\tC\t60\tPASS\t.\tGT\t1/1'],
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n')
        prepare(self.args);seq=self.consensus()
        self.assertEqual(len(seq),400) # short sequence excluded entirely
        self.assertEqual(seq[:10],'N'*10)
        self.assertEqual(seq[20],'R');self.assertEqual(seq[299],'C')
        self.assertEqual(seq[45:57],'N'*12)
        for lo,hi in [(100,110),(200,220),(390,400)]:self.assertEqual(seq[lo:hi],'N'*(hi-lo))
        self.assertEqual(seq[79],'N');self.assertEqual(seq[149],'N')
        self.assertEqual(seq[30],'A') # known homozygote inside confident BED
        qc=json.loads((self.root/'out/input_qc.json').read_text())
        self.assertEqual(qc['heterozygous_snps'],1)
    def test_read_absence_does_not_become_homozygous_reference(self):
        self.args.mode='hifi'
        d=self.root/'depth.tsv';d.write_text('#rname\tstartpos\tendpos\tmeandepth\nlong\t1\t400\t8\n')
        self.args.depth=str(d)
        fmt='GT:DP:AD:GQ:PL'
        self.vcf([f'21\t.\tA\tG\t60\tPASS\t.\t{fmt}\t0/1:8:4,4:60:50,0,50',
                  f'22\t.\tA\t<*>\t60\tPASS\t.\t{fmt}\t0/0:8:8,0:.:0,50,100',
                  f'23\t.\tA\t<*>\t60\tPASS\t.\t{fmt}\t0/0:8:8,0:.:0,5,60',
                  f'24\t.\tA\tG\t60\tPASS\t.\t{fmt}\t0/1:1:1,0:60:50,0,50',
                  f'25\t.\tA\tG\t60\tPASS\t.\t{fmt}\t0/1:8:7,1:60:50,0,50'],
                 ''.join(f'##FORMAT=<ID={name},Number={number},Type={kind},Description="test">\n' for name,number,kind in
                         [('GT','1','String'),('DP','1','Integer'),('AD','R','Integer'),('GQ','1','Integer'),('PL','G','Integer')]))
        prepare(self.args);seq=self.consensus()
        self.assertEqual(seq[20:25],'RANNN')
        self.assertEqual(seq[:20],'N'*20);self.assertEqual(seq[25:],'N'*375)
        qc=json.loads((self.root/'out/input_qc.json').read_text());self.assertEqual(qc['callable_bp'],2)
    def test_wrong_reference_hash_rejected(self):
        self.args.reference_sha256='f'*64
        with self.assertRaisesRegex(ValueError,'checksum'):prepare(self.args)
    def test_duplicate_variant_positions_are_masked_not_arbitrarily_chosen(self):
        conf=self.root/'confident.bed';conf.write_text('long\t0\t400\n');self.args.confident=str(conf)
        self.vcf(['21\t.\tA\tG\t60\tPASS\t.\tGT\t0/1','21\t.\tA\tC\t60\tPASS\t.\tGT\t0/1'],
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n')
        prepare(self.args)
        self.assertEqual(self.consensus()[15:26],'N'*11)
        self.assertEqual(json.loads((self.root/'out/input_qc.json').read_text())['duplicate_position_records_masked'],1)
    def test_poor_callability_explicit_not_fitted(self):
        out=Path(self.args.out);out.mkdir()
        (out/'input_qc.json').write_text(json.dumps(dict(callable_bp=0,callable_fraction=0,heterozygous_snps=0)))
        (out/'input.psmcfa').write_text('>c\nNNNN\n')
        (out/'contig_qc.tsv').write_text('contig\tcallable_bp\nc\t0\n')
        finalize(argparse.Namespace(out=str(out),min_callable_bp=10,min_callable_fraction=.1,min_hets=1,bin_size=100))
        qc=json.loads((out/'input_qc.json').read_text());self.assertEqual(qc['status'],'not_fitted')


if __name__=='__main__':unittest.main()
