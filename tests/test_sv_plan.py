import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from sv_plan import plan


class SVPlanTests(unittest.TestCase):
    def data(self):
        def assembly(name, tax='1'):
            return dict(id=name, taxid=tax, completion='finalized', chromosome_scale=True,
                        representation_eligible=True, reference_contigs=['chr1_1', 'chr2_1'], fasta=name + '.fa')
        return dict(assemblies=[assembly('reference'), assembly('query')],
                    cohorts=[dict(taxid='1', reference_id='reference', missing=[], status='disabled_but_eligible')])

    def test_graph_disabled_still_plans(self):
        result = plan(self.data())
        self.assertEqual(len(result['pairs']), 1)
        self.assertEqual(result['pairs'][0]['query'], 'query')
        self.assertEqual(result['pairs'][0]['chromosomes'], ['chr1_1', 'chr2_1'])

    def test_no_silent_intersection(self):
        data = self.data()
        data['assemblies'][1]['reference_contigs'] = ['chr1_1']
        result = plan(data)
        self.assertFalse(result['pairs'])
        self.assertEqual(result['outcomes'][0]['reason'], 'chromosome_sets_differ')

    def test_composites_and_fragmented_skip(self):
        for field, value, reason in [
                ('reference_contigs', ['chr1_1+chr2_1'], 'ambiguous_or_unharmonized_chromosome_names'),
                ('chromosome_scale', False, 'not_chromosome_scale')]:
            data = self.data()
            data['assemblies'][1][field] = value
            result = plan(data)
            self.assertEqual(result['outcomes'][0]['reason'], reason)

    def test_species_reference_cannot_cross(self):
        data = self.data()
        data['assemblies'][0]['taxid'] = '2'
        with self.assertRaises(ValueError):
            plan(data)

    def test_incomplete_species_withheld(self):
        data = self.data()
        data['cohorts'][0]['missing'] = ['absent']
        self.assertFalse(plan(data)['pairs'])


if __name__ == '__main__':
    unittest.main()
