import copy
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'py_scripts'))
from pangenome_input_audit import audit


class GraphInputAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        fa = Path(self.tmp.name)/'test.fa'
        fa.write_text('>chr1\nACGT\n')
        self.cohort = dict(taxid='1', min_individuals=2, reference_id='a', reference_name='REF_a', reference_contigs=['chr1'],
            members=[dict(id='a', sample='a', graph_name='REF_a', taxid='1', staged_fasta=str(fa)),
                     dict(id='b', sample='b', graph_name='IND_b', taxid='1', staged_fasta=str(fa))])

    def tearDown(self):
        self.tmp.cleanup()

    def test_identity_and_checksum(self):
        rows = audit(self.cohort)
        self.assertEqual(rows[0]['bases'], 4)
        self.assertTrue(rows[0]['is_reference'])
        self.assertEqual(rows[0]['sha256'], rows[1]['sha256'])

    def test_reference_contig_absent(self):
        self.cohort['reference_contigs'] = ['absent']
        with self.assertRaisesRegex(ValueError, 'absent from FASTA'):
            audit(self.cohort)

    def test_cross_species(self):
        self.cohort['members'][1]['taxid'] = '2'
        with self.assertRaisesRegex(ValueError, 'Cross-species'):
            audit(self.cohort)

    def test_haplotypes_do_not_create_individuals(self):
        self.cohort['members'][1]['sample'] = 'a'
        with self.assertRaisesRegex(ValueError, 'Insufficient biological'):
            audit(self.cohort)

    def test_reference_graph_name_mismatch(self):
        self.cohort['reference_name'] = 'wrong'
        with self.assertRaisesRegex(ValueError, 'identity disagree'):
            audit(self.cohort)

    def test_duplicate_sequence_names(self):
        Path(self.cohort['members'][0]['staged_fasta']).write_text('>chr1\nACGT\n>chr1\nA\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate FASTA'):
            audit(self.cohort)


if __name__ == '__main__':
    unittest.main()
