"""Native equivalence of whole-scaffold parallel and serial HiFi calling."""
import importlib.util
import json
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('read_call',Path(__file__).resolve().parents[1]/'py_scripts/psmc_read_call.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
try:import pysam
except ImportError:pysam=None


class ParallelReadCall(unittest.TestCase):
    @unittest.skipUnless(pysam and shutil.which('bcftools') and shutil.which('samtools'),'Requires native calling environment')
    def test_parallel_calls_equal_serial_in_fai_order_including_empty_scaffold(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rng=random.Random(37);names=['chrZ','chrA','empty','short']
            seqs={n:''.join(rng.choices('ACGT',k=1200 if n!='short' else 50)) for n in names}
            fa=root/'ref.fa';fa.write_text(''.join('>'+n+'\n'+s+'\n' for n,s in seqs.items()));pysam.faidx(str(fa))
            bam=root/'reads.bam'
            with pysam.AlignmentFile(str(bam),'wb',header=dict(HD=dict(VN='1.6',SO='coordinate'),
                SQ=[dict(SN=n,LN=len(s)) for n,s in seqs.items()],RG=[dict(ID='reads',SM='individual')])) as dest:
                for tid,name in enumerate(names[:2]):
                    for start in range(0,1000,25):
                        for rep in range(12):
                            seq=list(seqs[name][start:start+150])
                            if start<=500<start+150 and rep%2:seq[500-start]=next(b for b in 'ACGT' if b!=seq[500-start])
                            read=pysam.AlignedSegment();read.query_name=f'{name}_{start}_{rep}';read.query_sequence=''.join(seq)
                            read.flag=0;read.reference_id=tid;read.reference_start=start;read.mapping_quality=60
                            read.cigar=[(0,150)];read.query_qualities=pysam.qualitystring_to_array('I'*150);read.set_tag('RG','reads');dest.write(read)
            pysam.index(str(bam));bed=root/'all.bed';bed.write_text(''.join(f'{n}\t0\t1200\n' for n in names[:3]))
            serial=root/'serial.vcf.gz'
            with (root/'serial.log').open('w') as log:
                producer=subprocess.Popen(['bcftools','mpileup','-Ou','-f',str(fa),'-q','20','-Q','20','-R',str(bed),'-d','10000','-a','FORMAT/DP,FORMAT/AD',str(bam)],stdout=subprocess.PIPE,stderr=log)
                subprocess.run(['bcftools','call','--ploidy','2','-m','-A','--keep-unseen-allele','-a','GQ','-Oz','-o',str(serial)],stdin=producer.stdout,stderr=log,check=True)
                producer.stdout.close();self.assertEqual(producer.wait(),0)
            subprocess.run(['bcftools','index','--tbi',str(serial)],check=True)
            from argparse import Namespace
            out=root/'parallel';mod.run(Namespace(reference=str(fa),bam=str(bam),out=str(out),min_contig_bp=1000,mapq=20,baseq=20,max_pileup_depth=10000,cpus=4,workers=8))
            with pysam.VariantFile(str(serial)) as a,pysam.VariantFile(str(out/'diploid.vcf.gz')) as b:
                left=[str(r) for r in a];right=[str(r) for r in b]
            self.assertTrue(left);self.assertEqual(left,right)
            q=json.loads((out/'calling_manifest.json').read_text());self.assertEqual(q['workers'],2)
            self.assertEqual([r['scaffold'] for r in q['scaffolds']],names[:3])
            self.assertTrue((out/'diploid.vcf.gz.tbi').is_file())


if __name__=='__main__':unittest.main()
