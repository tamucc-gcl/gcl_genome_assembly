import gzip
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'py_scripts'))
from pangenome_output_audit import audit, REQUIRED


class OutputAuditTests(unittest.TestCase):
    def test_missing_required_has_role_specific_error(self):
        rows, errors = audit({})
        self.assertEqual(len(errors), len(REQUIRED))
        self.assertTrue(any('variant_gbz' in e for e in errors))

    def test_header_only_vcf_is_valid_empty_callset(self):
        with tempfile.TemporaryDirectory() as directory:
            roles = {}
            for role in REQUIRED:
                path = Path(directory)/role
                path.write_bytes(b'placeholder')
                roles[role] = str(path)
            with gzip.open(roles['variants_standard'], 'wt') as handle:
                handle.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            rows, errors = audit(roles)
            self.assertEqual(errors, [])
            self.assertEqual(next(r for r in rows if r['role']=='variants_standard')['graph'], 'gref_clip')
            Path(roles['variants_standard']).write_bytes(b'not gzip')
            self.assertTrue(any('variants_standard: invalid' in e for e in audit(roles)[1]))

    def test_optional_absence_does_not_add_required_error(self):
        rows, errors = audit({'variants_raw': ''})
        self.assertEqual(len(errors), len(REQUIRED))
        self.assertFalse(next(r for r in rows if r['role']=='variants_raw')['required'])


if __name__ == '__main__':
    unittest.main()
