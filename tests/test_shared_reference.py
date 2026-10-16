import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from shared_reference import verify


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.fa = Path(self.temp.name) / 'final.fa'
        self.fai = Path(self.temp.name) / 'final.fa.fai'
        self.fa.write_bytes(b'>chr1\nACGT\n>small\nNN\n')
        self.fai.write_text('chr1\t4\t6\t4\t5\nsmall\t2\t18\t2\t3\n')
        self.ref = dict(id='one_hap1', taxid='1', reference_contigs=['chr1'])

    def test_checksum_covers_all_final_records(self):
        result = verify(self.ref, self.fa, self.fai)
        self.assertEqual(result['sha256'], hashlib.sha256(self.fa.read_bytes()).hexdigest())
        self.assertEqual(result['bases'], 6)
        self.assertEqual(result['sequence_count'], 2)
        self.assertEqual(result['coordinate_stage'], 'final')

    def test_stale_index_rejected(self):
        self.fai.write_text('chr1\t3\t6\t4\t5\nsmall\t2\t18\t2\t3\n')
        with self.assertRaisesRegex(ValueError, 'disagree'):
            verify(self.ref, self.fa, self.fai)

    def test_changed_bases_change_identity_even_with_same_dictionary(self):
        before = verify(self.ref, self.fa, self.fai)
        self.fa.write_bytes(self.fa.read_bytes().replace(b'ACGT', b'AAAA'))
        after = verify(self.ref, self.fa, self.fai)
        self.assertNotEqual(before['sha256'], after['sha256'])
        self.assertEqual(before['sequences'], after['sequences'])

    def test_duplicate_index_names_rejected(self):
        self.fai.write_text('chr1\t4\nchr1\t4\n')
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            verify(self.ref, self.fa, self.fai)


if __name__ == '__main__':
    unittest.main()
