#!/usr/bin/env python3
"""Configure and verify the staged cohort acceptance run on existing Crest inputs."""
import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'py_scripts'))
from contamination_review import digest


def table(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def configure(project, preflight, out, expected_count=10):
    metadata = json.loads((preflight / 'preflight.json').read_text())
    inputs = table(preflight / 'inputs.tsv')
    if (metadata['status'] != 'PASS' or metadata['totals']['assemblies'] != expected_count or
            len(inputs) != expected_count or len({row['assembly'] for row in inputs}) != expected_count or
            Path(metadata['project']).resolve() != project.resolve()):
        raise ValueError('Expected a successful whole-cohort preflight for this project')
    expected = {}
    for row in inputs:
        assembly = row['assembly']
        source_path = preflight / 'assemblies' / assembly / 'queries/source.json'
        manifest = json.loads((preflight / 'assemblies' / assembly / 'preflight_manifest.json').read_text())
        if digest(source_path) != manifest['outputs']['queries/source.json']:
            raise ValueError('Changed query manifest: ' + assembly)
        source = json.loads(source_path.read_text())
        if digest(Path(row['fasta'])) != source['source_sha256']:
            raise ValueError('Changed source FASTA: ' + assembly)
        expected[assembly] = {key: source[key] for key in ('source_sha256', 'source_bp', 'sequences')}
    protein = project / 'comparisons/contamination-full-protein-01/db/diamond/uniprot_reference_proteomes/2026_04'
    taxonomy = project / 'comparisons/contamination-full-protein-01/db/ncbi_taxonomy/2026_04'
    nt = project / 'comparisons/contamination-ctlk-whole-02/db/blast/nt/2026-11-08/database/nt'
    for path in (project / 'data/assembly_samplesheet.csv', project / 'data/hic_readsets.csv',
                 protein / 'reference_proteomes.dmnd', protein / 'protein_database.json',
                 taxonomy / 'nodes.dmp', taxonomy / 'names.dmp', taxonomy / 'merged.dmp', Path(str(nt) + '.nal')):
        if not path.is_file() or not path.stat().st_size:
            raise ValueError('Missing retained input: ' + str(path))
    with (project / 'data/assembly_samplesheet.csv').open(newline='') as handle:
        samples = [row['sample_id'] for row in csv.DictReader(handle)]
    if (len(samples) != len(set(samples)) or
            set(samples) != {assembly.rsplit('_hap', 1)[0] for assembly in expected}):
        raise ValueError('Sample sheet membership differs from the preflight cohort')
    params = dict(sample_sheet=str(project / 'data/assembly_samplesheet.csv'),
                  hic_readsets=str(project / 'data/hic_readsets.csv'), publish_dir_mode='link',
                  taxid=373251, ploidy=2, qc_mode='none', run_purge_dups=True,
                  inspector_contig_min_contig_bp=25000, inspector_contig_assemblyerror_min_bp=200000,
                  inspector_contig_min_assembly_error_bp=5000, inspector_contig_max_assembly_error_bp=50000000,
                  run_decon_contigs=True, run_contamination_evidence=True, contamination_search_mode='staged',
                  contamination_read_support=True, contamination_review=None, stop_after_contig_review=True,
                  contamination_protein_db_dir=str(protein), contamination_protein_taxonomy_dir=str(taxonomy),
                  contamination_nucleotide_db=str(nt), contamination_targeted_blastn=False,
                  contamination_blastx_max_forks=56, contamination_megablast_execution='nodes',
                  contamination_megablast_nodes=4, contamination_megablast_slots=8,
                  contamination_megablast_threads=8, contamination_megablast_time='48h',
                  contamination_megablast_scratch='/scratch')
    for key, value in metadata['settings'].items():
        if key != 'chunk_bp':
            params['contamination_sample_' + key] = value
        elif value != 5000000:
            raise ValueError('Unexpected protein chunk size')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'params.json').write_text(json.dumps(params, indent=2) + '\n')
    (out / 'expected.json').write_text(json.dumps(expected, indent=2) + '\n')
    print(f'Configured {len(expected)} assemblies; retained BAMs enabled; no removal sheet; stop before scaffolding.', flush=True)


def verify(expected_path, results, out):
    expected = json.loads(expected_path.read_text())
    directories = {path.name: path for path in (results / 'assembly/contamination').iterdir()
                   if path.is_dir() and (path / 'review/evidence.json').is_file()}
    if directories.keys() != expected.keys():
        raise ValueError('Report assembly set differs from expected cohort')
    summaries = []
    for assembly, directory in sorted(directories.items()):
        review = directory / 'review'
        metadata = json.loads((review / 'evidence.json').read_text())
        if metadata['assembly'] != assembly:
            raise ValueError('Report identity mismatch: ' + assembly)
        for key, value in expected[assembly].items():
            if metadata[key] != value:
                raise ValueError(f'{assembly} source differs from preflight: {key}')
        if metadata['search_mode'] != 'staged' or metadata['coverage'] != 'provided-bam':
            raise ValueError('Expected staged evidence with BAM-derived coverage: ' + assembly)
        plots = table(review / 'locus_plot_hits.tsv')
        ids = {row['action_id'] for row in plots}
        depth = table(review / 'locus_read_depth.tsv')
        boundaries = table(review / 'locus_read_boundaries.tsv')
        if {row['action_id'] for row in depth} != ids:
            raise ValueError('Missing BAM depth evidence: ' + assembly)
        if len(boundaries) != 2 * len(ids) or {row['action_id'] for row in boundaries} != ids:
            raise ValueError('Missing boundary read evidence: ' + assembly)
        for identity in ids:
            if {row['edge'] for row in boundaries if row['action_id'] == identity} != {'Left', 'Right'}:
                raise ValueError('Incomplete boundary pair: ' + identity)
            figure = directory / 'diagnostics/loci' / (identity + '.png')
            if not figure.is_file() or not figure.stat().st_size:
                raise ValueError('Missing locus figure: ' + identity)
        for name in ('blob.png', 'cumulative.png'):
            figure = directory / 'diagnostics' / name
            if not figure.is_file() or not figure.stat().st_size:
                raise ValueError('Missing diagnostics: ' + assembly)
        if any(row['selected'] != 'NO' for row in table(review / 'review.tsv')):
            raise ValueError('Unexpected selected removal in validation: ' + assembly)
        summaries.append(dict(assembly=assembly, source_bp=metadata['source_bp'],
                              protein_query_bp=metadata['protein_searched_bp'],
                              nucleotide_query_bp=metadata['megablast_searched_bp'],
                              protein_chunks=metadata['chunks'], locus_panels=len(ids)))
    completion = json.loads((results / 'pipeline/contamination_megablast_completion.json').read_text())
    if completion['missing']:
        raise ValueError('Nucleotide chunks are incomplete')
    out.write_text(json.dumps(dict(status='PASS', assemblies=summaries,
                                  skipped_nucleotide_workers=completion['skipped_workers'],
                                  note='Filtered BAM coverage may select different windows than the table-only preflight.'), indent=2) + '\n')
    print(f'PASS: {len(summaries)} staged reports, BAM read support, figures and complete nucleotide chunks', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    c = commands.add_parser('configure')
    for name in ('project', 'preflight', 'out'):
        c.add_argument('--' + name, type=Path, required=True)
    v = commands.add_parser('verify')
    for name in ('expected', 'results', 'out'):
        v.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'configure':
        configure(args.project.resolve(strict=True), args.preflight.resolve(strict=True), args.out)
    else:
        verify(args.expected, args.results, args.out)


if __name__ == '__main__':
    main()
