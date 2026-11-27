import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_contamination_staged_preflight as fixtures
from prepare_contamination_staged_cohort import configure, verify, validate_sample_membership
import contamination_staged_preflight as preflight


class CohortTests(fixtures.PreflightTests):
    test_complete_census_reuses_identical_queries_and_excludes_fastas_from_archive = None
    test_rejects_incomplete_cohort_mismatched_coverage_and_changed_queries = None
    test_two_process_cli_needs_no_alignment_executables = None

    def test_configuration_and_bam_acceptance(self):
        preflight.run(self.args)
        for relative in ('data/assembly_samplesheet.csv', 'data/hic_readsets.csv',
                         'comparisons/contamination-full-protein-01/db/diamond/uniprot_reference_proteomes/2026_04/reference_proteomes.dmnd',
                         'comparisons/contamination-full-protein-01/db/diamond/uniprot_reference_proteomes/2026_04/protein_database.json',
                         'comparisons/contamination-full-protein-01/db/ncbi_taxonomy/2026_04/nodes.dmp',
                         'comparisons/contamination-full-protein-01/db/ncbi_taxonomy/2026_04/names.dmp',
                         'comparisons/contamination-full-protein-01/db/ncbi_taxonomy/2026_04/merged.dmp',
                         'comparisons/contamination-ctlk-whole-02/db/blast/nt/2026-11-08/database/nt.nal'):
            path = self.project / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('fixture')
        (self.project / 'data/assembly_samplesheet.csv').write_text(
            'sample_id,hifi_bam,hic_r1,hic_r2\nfixture,fixture.bam,,\nhic_only,,R1.fq,R2.fq\n')
        run = self.project / 'cohort'
        configure(self.project, self.out, run, expected_count=2)
        params = json.loads((run / 'params.json').read_text())
        self.assertTrue(params['run_decon_contigs'] and params['contamination_read_support'])
        self.assertIsNone(params['contamination_review'])
        self.assertTrue(params['stop_after_contig_review'])
        self.assertEqual((params['contamination_blastx_max_forks'], params['contamination_megablast_nodes']), (56, 4))
        expected = json.loads((run / 'expected.json').read_text())
        results = run / 'results'
        for assembly, source in expected.items():
            directory = results / 'assembly/contig/decontam/blobtools' / assembly
            review = directory / 'review'
            review.mkdir(parents=True)
            metadata = dict(source, assembly=assembly, search_mode='staged', coverage='provided-bam',
                            protein_searched_bp=1000, megablast_searched_bp=500, chunks=1)
            (review / 'evidence.json').write_text(json.dumps(metadata))
            (review / 'locus_plot_hits.tsv').write_text('action_id\nA1\n')
            (review / 'locus_read_depth.tsv').write_text('action_id\nA1\n')
            (review / 'locus_read_boundaries.tsv').write_text('action_id\tedge\nA1\tLeft\nA1\tRight\n')
            (review / 'review.tsv').write_text('selected\nNO\n')
            figures = directory / 'diagnostics/loci'
            figures.mkdir(parents=True)
            for path in (figures / 'A1.png', figures.parent / 'blob.png', figures.parent / 'cumulative.png'):
                path.write_bytes(b'figure transport fixture')
        (results / 'pipeline').mkdir()
        (results / 'pipeline/contamination_megablast_completion.json').write_text(json.dumps(dict(missing=[], skipped_workers=[])))
        cohort = results / 'assembly/contig/decontam/report'
        cohort.mkdir()
        (cohort / 'report.json').write_text(json.dumps(dict(fcs_enabled=True, blobtools_enabled=True,
            assemblies=[dict(assembly=name) for name in expected])))
        for name in ('cohort_gc_coverage.png', 'cohort_gc_coverage.pdf', 'cohort_gc_coverage_scope.tsv'):
            (cohort / name).write_text('figure transport fixture')
        fcs = results / 'assembly/contig/decontam/fcs'
        fcs.mkdir()
        for assembly in expected:
            (fcs / f'{assembly}.fcs_gx_report.txt').write_text('#seq_id\tstart_pos\tend_pos\taction\n')
        verify(run / 'expected.json', results, run / 'acceptance.json')
        self.assertEqual(json.loads((run / 'acceptance.json').read_text())['status'], 'PASS')
        (results / 'assembly/contig/decontam/blobtools/fixture_hap1/review/locus_read_boundaries.tsv').write_text('action_id\tedge\nA1\tLeft\n')
        with self.assertRaisesRegex(ValueError, 'Missing boundary'):
            verify(run / 'expected.json', results, run / 'invalid.json')

    def test_membership_keeps_hic_only_rows_but_rejects_assembly_mismatches(self):
        sheet = self.project / 'samples.csv'
        sheet.write_text('\ufeffsample_id, hifi_bam, sr_r1, sr_r2, hic_r1, hic_r2\n'
                         'hifi,reads.bam,,,,\nshotgun,,R1.fq,R2.fq,,\nhic_only,,,,H1.fq,H2.fq\n\n', encoding='utf-8')
        self.assertEqual(validate_sample_membership(sheet, ['hifi_hap1', 'hifi_hap2', 'shotgun_hap1']), ['hic_only'])
        with self.assertRaisesRegex(ValueError, "unexpected=\\['shotgun'\\]"):
            validate_sample_membership(sheet, ['hifi_hap1', 'hifi_hap2'])
        with self.assertRaisesRegex(ValueError, "missing=\\['absent'\\]"):
            validate_sample_membership(sheet, ['hifi_hap1', 'shotgun_hap1', 'absent_hap1'])
        sheet.write_text('sample_id,hifi_bam\nhifi,a.bam\nhifi,b.bam\n')
        with self.assertRaisesRegex(ValueError, 'duplicate sample_id'):
            validate_sample_membership(sheet, ['hifi_hap1'])


if __name__ == '__main__':
    unittest.main()
