"""Known-answer native Dipcall length tests; mask growth must have matching calls."""
import gzip
import hashlib
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

try:import pysam
except ImportError:pysam=None


class NativeAlignmentLengths(unittest.TestCase):
    @unittest.skipUnless(pysam and shutil.which('run-dipcall'),'Requires native Dipcall environment')
    def test_short_alignment_masks_and_variant_calls_are_enabled_together(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rng=random.Random(613);ref=''.join(rng.choices('ACGT',k=250000))
            alternate=list(ref)
            positions=[1000,31000,50000,150000]
            for pos in positions:alternate[pos]=next(b for b in 'ACGT' if b!=ref[pos])
            alternate=''.join(alternate)
            fa=root/'ref.fa';fa.write_text('>chr1\n'+ref+'\n');h1=root/'h1.fa';h1.write_text(fa.read_text())
            h2=root/'h2.fa';pieces=[(0,30000),(30000,45000),(45000,105000),(105000,250000)]
            h2.write_text(''.join(f'>piece{i}\n{alternate[lo:hi]}\n' for i,(lo,hi) in enumerate(pieces)))
            tool=Path(shutil.which('run-dipcall')).parent
            fingerprints={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [tool/'paftools.js',tool/'dipcall-aux.js']}
            expected={50000:[False,False,True,True],20000:[True,False,True,True],10000:[True,True,True,True]}
            for minimum,enabled in expected.items():
                work=root/str(minimum);work.mkdir();(work/'reference.fa').symlink_to(fa);pysam.faidx(str(work/'reference.fa'))
                subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'py_scripts/psmc_dipcall.py'),
                    '--preset','asm20','--min-alignment-bp',str(minimum),'--threads-per-hap','1','--haplotypes',str(h1),str(h2)],cwd=work,check=True)
                intervals=[line.split() for line in (work/'diploid.dip.bed').read_text().splitlines()]
                for pos,want in zip(positions,enabled):
                    masked=any(int(r[1])<=pos<int(r[2]) for r in intervals)
                    self.assertEqual(masked,want,(minimum,pos,intervals))
                subprocess.run(['bcftools','index','--tbi',str(work/'diploid.dip.vcf.gz')],check=True)
                with pysam.VariantFile(str(work/'diploid.dip.vcf.gz')) as vcf:
                    heterozygotes={r.start for r in vcf if r.samples[0]['GT'] in [(0,1),(1,0)] and not (set(r.filter)-{'PASS','.'})}
                for pos,want in zip(positions,enabled):self.assertEqual(pos in heterozygotes,want,(minimum,pos,heterozygotes))
            for name,value in fingerprints.items():self.assertEqual(hashlib.sha256((tool/name).read_bytes()).hexdigest(),value)


if __name__=='__main__':unittest.main()
