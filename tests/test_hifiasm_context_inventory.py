"""Bounded trace lookup safeguards; no assembly tools or analyses."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/comparisons'))
from inventory_hifiasm_context import trace_task_roots


class ContextInventoryTests(unittest.TestCase):
    def test_only_exact_sample_and_unique_task_are_resolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root/'b4/4c0468abc').mkdir(parents=True)
            trace = root/'trace.tsv'
            trace.write_text('name\thash\nCONTIG_ASSEMBLY:HIFIASM (fish)\tb4/4c0468\n'
                             'CONTIG_ASSEMBLY:HIFIASM (other)\tzz/unsafe\n')
            roots, audit = trace_task_roots(trace, root, 'fish')
            self.assertEqual(roots, [root/'b4/4c0468abc'])
            self.assertEqual(audit[0]['status'], 'resolved')
            (root/'b4/4c0468def').mkdir()
            self.assertEqual(trace_task_roots(trace, root, 'fish')[0], [])
            self.assertEqual(trace_task_roots(trace, root, 'fish')[1][0]['status'], 'ambiguous')

    def test_missing_task_and_invalid_hash_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); trace = root/'trace.tsv'
            trace.write_text('name\thash\nCONTIG_ASSEMBLY:HIFIASM (fish)\tb4/4c0468\n')
            self.assertEqual(trace_task_roots(trace, root, 'fish')[1][0]['status'], 'missing')
            trace.write_text('name\thash\nCONTIG_ASSEMBLY:HIFIASM (fish)\t../../outside\n')
            with self.assertRaises(ValueError): trace_task_roots(trace, root, 'fish')


if __name__ == '__main__':
    unittest.main()
