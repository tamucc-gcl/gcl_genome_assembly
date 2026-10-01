import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from pangenome_variant_audit import sample_mapping


class VariantIdentityTests(unittest.TestCase):
    def setUp(self):
        self.rows = [dict(id='a1', sample='a', graph_name='REF_a1'),
                     dict(id='a2', sample='a', graph_name='IND_a.2'),
                     dict(id='b1', sample='b', graph_name='IND_b.1'),
                     dict(id='b2', sample='b', graph_name='IND_b.2')]

    def test_reference_is_not_an_extra_individual(self):
        result = sample_mapping(self.rows, ['IND_b', 'REF_a1', 'IND_a'])
        self.assertEqual(len(result), 3)
        self.assertEqual(len({r[1] for r in result}), 2)
        self.assertEqual(result[0][3], 2)

    def test_missing_and_unexpected_columns_fail(self):
        for samples in (['IND_b', 'IND_a'], ['IND_b', 'REF_a1', 'IND_a', 'synthetic'],
                        ['IND_b', 'REF_a1', 'IND_a', 'IND_a']):
            with self.assertRaises(ValueError):
                sample_mapping(self.rows, samples)


if __name__ == '__main__':
    unittest.main()
