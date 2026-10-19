import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'py_scripts'))
from pangenome_regional_context import prepare


class RegionalContextTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.sharing = self.root/'sharing.tsv'
        self.sharing.write_text('scope\thaplotype\tchromosome\n'
            'haplotype\tA\t\nhaplotype\tB\t\n'
            'chromosome_haplotype\tA\tchr4\n'
            'chromosome_haplotype\tA\tcomposite:chr4_1+chr14_1\n'
            'chromosome_haplotype\tB\tcomposite:chr4_1+chr14_1\n')
        self.eligibility = dict(assemblies=[
            dict(id='A', taxid='1', chromosome_scale=True, reference_contigs=['chr4_1+chr14_1','chr4_2']),
            dict(id='B', taxid='1', chromosome_scale=False, reference_contigs=['chr4_1+chr14_1'])])

    def test_scope_comes_from_final_eligibility_not_composite_spelling(self):
        groups, presence = prepare(self.sharing, self.eligibility, '1', self.root/'out')
        comps = {r['haplotype']:r for r in groups if r['kind']=='composite'}
        self.assertEqual(comps['A']['chromosome_scale_composite'], 'yes')
        self.assertEqual(comps['B']['chromosome_scale_composite'], 'no')
        keyed = {(r['haplotype'],r['chromosome']):r for r in presence}
        self.assertEqual(keyed['A','chr4']['status'], 'separate_and_composite')
        self.assertEqual(keyed['A','chr14']['status'], 'composite_only')
        self.assertEqual(keyed['B','chr14']['status'], 'other_composite_only')

    def test_other_species_cannot_supply_context(self):
        with self.assertRaisesRegex(ValueError, 'species eligibility'):
            prepare(self.sharing, self.eligibility, '2', self.root/'out')

    def test_small_scaffold_in_large_assembly_is_still_supplementary(self):
        self.eligibility['assemblies'][0]['reference_contigs'] = ['chr4_2']
        groups, _ = prepare(self.sharing, self.eligibility, '1', self.root/'out')
        self.assertTrue(all(r['chromosome_scale_composite']=='no' for r in groups))

    def test_source_composite_absent_from_graph_is_not_hidden(self):
        self.sharing.write_text(self.sharing.read_text().replace('chromosome_haplotype\tA\tcomposite:chr4_1+chr14_1\n',''))
        groups, presence = prepare(self.sharing, self.eligibility, '1', self.root/'out')
        missing=next(r for r in groups if r['haplotype']=='A' and r['kind']=='composite')
        self.assertEqual(missing['graph_group_present'],'no')
        keyed={(r['haplotype'],r['chromosome']):r for r in presence}
        self.assertEqual(keyed['A','chr4']['status'],'separate_and_composite_missing')
        self.assertEqual(keyed['A','chr14']['status'],'composite_missing_from_graph')


if __name__ == '__main__':
    unittest.main()
