"""Cleanup acceptance preserves base identity while accounting for GX case conversion."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'py_scripts'))
spec=importlib.util.spec_from_file_location('adaptor_acceptance',ROOT/'tests/fcs_adaptor_acceptance.py')
acceptance=importlib.util.module_from_spec(spec)
spec.loader.exec_module(acceptance)


class CleanupAcceptance(unittest.TestCase):
    def setUp(self):
        self.source={'unchanged':'aCgTnN','seq_00018':'a'*522}
        self.cleaned={'unchanged':'ACGTNN'}
        self.actions=[dict(accession='seq_00018',action='ACTION_EXCLUDE')]
        for i in range(1,11):
            name=f'seq_{i:05d}'
            sequence='acGt'*25100
            self.source[name]=sequence
            lo,hi=(1,58) if i<=5 else (100001,100058)
            self.actions.append(dict(accession=name,action='ACTION_TRIM',range=f'{lo}..{hi}'))
            if lo>1:self.cleaned[f'{name}~1..{lo-1}']=sequence[:lo-1].upper()
            self.cleaned[f'{name}~{hi+1}..{len(sequence)}']=sequence[hi:].upper()

    def check(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'source.fa'
            source.write_text('fixture')
            out=Path(temp)/'acceptance.json'
            with patch.object(acceptance,'fasta',side_effect=[self.source,self.cleaned]), \
                 patch.object(acceptance,'adaptor_rows',side_effect=[self.actions,[]]):
                acceptance.check(source,'initial','clean','repeat',out)
            return json.loads(out.read_text())

    def test_uppercase_output_including_unchanged_records(self):
        result=self.check()
        self.assertEqual(result['status'],'PASS')
        self.assertEqual(result['removed_bp'],1102)
        self.assertEqual(result['sequence_case'],'uppercase')

    def test_base_substitution_still_fails(self):
        self.cleaned['unchanged']='TCGTNN'
        with self.assertRaisesRegex(ValueError,"sequence_mismatches=.*unchanged"):
            self.check()

    def test_wrong_trim_boundary_still_fails(self):
        name='seq_00001~59..100400'
        self.cleaned[name]=self.cleaned[name][1:]
        with self.assertRaisesRegex(ValueError,'sequence_mismatches=.*seq_00001'):
            self.check()


if __name__=='__main__':unittest.main()
