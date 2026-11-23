#!/usr/bin/env python3
"""Bind a full-haplotype staged test to its completed preparation-only census."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'py_scripts'))
from contamination_review import digest


def configure(project, preflight, assembly, out):
    metadata = json.loads((preflight / 'preflight.json').read_text())
    if metadata['status'] != 'PASS' or metadata['alignment_jobs_launched'] != 0:
        raise ValueError('A successful preparation-only preflight is required')
    directory = preflight / 'assemblies' / assembly
    manifest = json.loads((directory / 'preflight_manifest.json').read_text())
    inputs = manifest['inputs']
    if inputs['assembly'] != assembly or Path(metadata['project']).resolve() != project.resolve():
        raise ValueError('Preflight project or assembly differs from this test')
    for key in ('fasta', 'coverage'):
        if digest(Path(inputs[key])) != inputs[key + '_sha256']:
            raise ValueError('Changed preflight input: ' + inputs[key])
    source = directory / 'queries/source.json'
    if digest(source) != manifest['outputs']['queries/source.json']:
        raise ValueError('Changed preflight query manifest')
    expected = json.loads(source.read_text())
    if expected['source_sha256'] != inputs['fasta_sha256'] or expected['search_mode'] != 'staged':
        raise ValueError('Query manifest does not match its preflight source')
    settings = inputs['settings']
    if settings != metadata['settings'] or settings['chunk_bp'] != 5000000:
        raise ValueError('Unexpected or inconsistent preflight sampling settings')
    protein = project / 'comparisons/contamination-full-protein-01/db/diamond/uniprot_reference_proteomes/2026_04'
    taxonomy = project / 'comparisons/contamination-full-protein-01/db/ncbi_taxonomy/2026_04'
    nt = project / 'comparisons/contamination-ctlk-whole-02/db/blast/nt/2026-11-08/database/nt'
    for path in (protein / 'reference_proteomes.dmnd', protein / 'protein_database.json',
                 taxonomy / 'nodes.dmp', taxonomy / 'names.dmp', Path(str(nt) + '.nal')):
        if not path.is_file() or not path.stat().st_size:
            raise ValueError('Missing retained database file: ' + str(path))
    params = dict(fasta=inputs['fasta'], coverage_tsv=inputs['coverage'], assembly=assembly,
                  taxid=373251, stage='post_fcs', search_mode='staged', targeted_blastn=False,
                  protein_db_dir=str(protein), taxonomy_dir=str(taxonomy), nucleotide_db=str(nt),
                  contamination_megablast_execution='nodes', contamination_megablast_nodes=1,
                  contamination_megablast_slots=8, contamination_megablast_threads=8,
                  contamination_megablast_scratch='/scratch', contamination_megablast_time='24h',
                  contamination_blastx_max_forks=16)
    for key, value in settings.items():
        if key != 'chunk_bp':
            params['contamination_sample_' + key] = value
    # The sampler's short_contig_bp key maps directly to the assembly parameter.
    out.mkdir(parents=True, exist_ok=True)
    (out / 'params.json').write_text(json.dumps(params, indent=2) + '\n')
    (out / 'expected.json').write_text(json.dumps(expected, indent=2) + '\n')
    print(f"{assembly}: whole source {expected['source_bp']:,} bp; "
          f"{expected['search_query_bp']:,} protein-query bp in {expected['chunks']} chunks", flush=True)


def verify(expected, evidence):
    expected = json.loads(expected.read_text())
    evidence = json.loads(evidence.read_text())
    for field in ('assembly', 'source_sha256', 'source_bp', 'sequences', 'search_mode',
                  'search_queries', 'search_query_bp', 'chunks', 'sampling'):
        if expected[field] != evidence[field]:
            raise ValueError('Full-haplotype result differs from preflight: ' + field)
    print('PASS: full-haplotype source, sampling and query counts match the preflight', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    configure_parser = sub.add_parser('configure')
    configure_parser.add_argument('--project', type=Path, required=True)
    configure_parser.add_argument('--preflight', type=Path, required=True)
    configure_parser.add_argument('--assembly', choices=['Sde-CTlk_104_hap1'], default='Sde-CTlk_104_hap1')
    configure_parser.add_argument('--out', type=Path, required=True)
    verify_parser = sub.add_parser('verify')
    verify_parser.add_argument('--expected', type=Path, required=True)
    verify_parser.add_argument('--evidence', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'configure':
        configure(args.project.resolve(strict=True), args.preflight.resolve(strict=True), args.assembly, args.out)
    else:
        verify(args.expected, args.evidence)


if __name__ == '__main__':
    main()
