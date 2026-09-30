import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('sequence_context', Path(__file__).resolve().parents[1]/'py_scripts/chimera_sequence_context.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SequenceContextTests(unittest.TestCase):
    def test_unresolved_chromosome_scope_is_excluded(self):
        row = dict(chromosome_member='unresolved', candidate_verdict='REVIEW', evidence_only='yes')
        self.assertEqual(module.eligible([row]), [])

    def test_diagnostic_interval_is_included_without_cut_permission(self):
        row = dict(chromosome_member='yes', candidate_verdict='REVIEW', evidence_only='yes', callable='no')
        self.assertEqual(module.eligible([row]), [row])
        self.assertEqual(row['callable'], 'no')

    def test_midpoint_only_alignment_is_not_whole_interval_support(self):
        self.assertFalse(module.spans_interval(['q', '500000', '230000', '270000'], 160000, 330000))

    def test_synthetic_sister_alignment_brackets_interval(self):
        self.assertTrue(module.spans_interval(['q', '500000', '158981', '500000'], 160739, 330451))

    def test_flanks_are_required(self):
        self.assertFalse(module.spans_interval(['q', '500000', '160739', '330451'], 160739, 330451))


if __name__ == '__main__':
    unittest.main()
