import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_contamination_staged_preflight as fixtures
from prepare_contamination_staged_haplotype import configure, verify
import contamination_staged_preflight as preflight


class HaplotypeTest(fixtures.PreflightTests):
    # Use just the fixture setup; the preflight's own tests remain in their module.
    test_complete_census_reuses_identical_queries_and_excludes_fastas_from_archive = None
    test_rejects_incomplete_cohort_mismatched_coverage_and_changed_queries = None
    test_two_process_cli_needs_no_alignment_executables = None

    def test_bound_configuration_and_result_verification(self):
        preflight.run(self.args)
        for path in ('comparisons/contamination-full-protein-01/db/diamond/uniprot_reference_proteomes/2026_04/reference_proteomes.dmnd',
                     'comparisons/contamination-full-protein-01/db/diamond/uniprot_reference_proteomes/2026_04/protein_database.json',
                     'comparisons/contamination-full-protein-01/db/ncbi_taxonomy/2026_04/nodes.dmp',
                     'comparisons/contamination-full-protein-01/db/ncbi_taxonomy/2026_04/names.dmp',
                     'comparisons/contamination-ctlk-whole-02/db/blast/nt/2026-11-08/database/nt.nal'):
            target = self.project / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text('fixture')
        run = self.project / 'haplotype'
        configure(self.project, self.out, 'fixture_hap1', run)
        params = json.loads((run / 'params.json').read_text())
        self.assertEqual(params['contamination_megablast_execution'], 'nodes')
        self.assertEqual(params['contamination_megablast_nodes'], 1)
        self.assertEqual(params['contamination_sample_short_contig_bp'], 100000)
        expected = json.loads((run / 'expected.json').read_text())
        evidence = run / 'evidence.json'
        evidence.write_text(json.dumps(expected))
        verify(run / 'expected.json', evidence)
        expected['search_query_bp'] += 1
        evidence.write_text(json.dumps(expected))
        with self.assertRaisesRegex(ValueError, 'search_query_bp'):
            verify(run / 'expected.json', evidence)
        source = Path(params['fasta'])
        source.write_text(source.read_text() + 'A\n')
        with self.assertRaisesRegex(ValueError, 'Changed preflight input'):
            configure(self.project, self.out, 'fixture_hap1', run)


if __name__ == '__main__':
    unittest.main()
