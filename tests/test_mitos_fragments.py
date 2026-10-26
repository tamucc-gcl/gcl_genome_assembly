"""Preserve all fragmented organelle records and their matching MITOS annotations."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from Bio import SeqIO


class MitosFragmentTests(unittest.TestCase):
    def test_every_fragment_has_its_own_sequence_and_annotation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            fa=root/'input.fa';fa.write_text('>one\nACGTACGTACGT\n>two\nTTTTCCCC\n')
            for i,gene in enumerate(['cox1','trnF']):
                bed=root/'mitos'/str(i)/'result.bed';bed.parent.mkdir(parents=True)
                bed.write_text(f'record\t0\t3\t{gene}\t0\t+\n')
            command=[sys.executable,str(Path(__file__).resolve().parents[1]/'py_scripts/mitos_to_genbank.py'),
                     '--fasta',str(fa),'--bed-dir',str(root/'mitos'),'--sample_id','sample',
                     '--output',str(root/'out.gb'),'--combined-bed',str(root/'out.bed')]
            subprocess.run(command,check=True,capture_output=True)
            with (root/'out.gb').open() as handle:
                records=list(SeqIO.parse(handle,'genbank'))
            self.assertEqual([str(r.seq) for r in records],['ACGTACGTACGT','TTTTCCCC'])
            self.assertEqual([r.features[1].qualifiers['gene'][0] for r in records],['cox1','trnF'])
            self.assertTrue(all(r.annotations['topology']=='linear' for r in records))
            self.assertEqual([line.split()[0] for line in (root/'out.bed').read_text().splitlines()],['one','two'])
            (root/'mitos/1/result.bed').unlink()
            failed=subprocess.run(command,capture_output=True,text=True)
            self.assertNotEqual(failed.returncode,0)
            self.assertIn('annotation is incomplete',failed.stderr)


if __name__=='__main__':unittest.main()
