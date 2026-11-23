#!/usr/bin/env python3
"""Prepare a whole-cohort staged-search census without databases or alignments."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import json
from pathlib import Path
import re
import sys
import tarfile
import tempfile
import time
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'py_scripts'))
from contamination_review import digest, prepare, write_table
from contamination_staged import sample
from contamination_evidence import table
from prepare_contamination_real_pilot import traced_task, subset_coverage


def discover(project, trace, expected):
    with trace.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    names = {}
    for row in rows:
        match = re.search(r'(?:^|:)FCS_CLEAN_GENOME \(([^)]+)\)$', row['name'])
        if match and row['status'] in {'CACHED', 'COMPLETED'} and row['exit'] == '0':
            names.setdefault(match[1], set()).add(row['name'])
    if len(names) != expected:
        raise ValueError(f'Expected {expected} post-FCS assemblies in {trace}; found {len(names)}: {sorted(names)}')
    result = []
    for assembly, tasks in sorted(names.items()):
        if not re.fullmatch('[A-Za-z0-9_.-]+', assembly) or assembly in {'.', '..'} or len(tasks) != 1:
            raise ValueError('Unsafe or ambiguous assembly identity: ' + assembly)
        fasta_task = traced_task(trace, project / 'work', next(iter(tasks)))
        fasta = (fasta_task / f'{assembly}.decontaminated.fasta').resolve(strict=True)
        candidates = {row['name'] for row in rows
                      if row['name'].endswith(f':ASSEMBLY_QC_CONTIG_DECONTAM:MAPPING_QC ({assembly})')
                      and row['status'] in {'CACHED', 'COMPLETED'} and row['exit'] == '0'}
        if len(candidates) != 1:
            raise ValueError(f'Expected one post-FCS QC coverage task for {assembly}; found {candidates}')
        coverage_task = traced_task(trace, project / 'work', candidates.pop())
        coverage = (coverage_task / f'{assembly}_mapping_stats/coverage.txt').resolve(strict=True)
        result.append(dict(assembly=assembly, fasta=str(fasta), coverage=str(coverage)))
    return result


def census(directory):
    source = json.loads((directory / 'queries/source.json').read_text())
    windows = table(directory / 'queries/windows.tsv')
    chunks = table(directory / 'queries/chunks.tsv')
    selected = [row for row in windows if row['status'] == 'search']
    bp = sum(int(row['end0']) - int(row['start0']) for row in selected)
    if (len(windows) != source['windows'] or len(selected) != source['search_queries'] or
            bp != source['search_query_bp'] or len(chunks) != source['chunks'] or
            sum(int(row['query_bp']) for row in chunks) != bp):
        raise ValueError('Prepared-query accounting mismatch: ' + str(directory))
    reasons = {}
    for row in selected:
        key = row['selection_reason']
        count, bases = reasons.get(key, (0, 0))
        reasons[key] = (count + 1, bases + int(row['end0']) - int(row['start0']))
    summary = dict(assembly=source['assembly'], source_bp=source['source_bp'], contigs=source['sequences'],
                   protein_query_bp=bp, protein_query_pct=100 * bp / source['source_bp'],
                   protein_windows=len(selected), protein_chunks=len(chunks), total_windows=len(windows),
                   genome_gc=source['sampling']['genome_gc'], depth_baseline=source['sampling']['depth_baseline'])
    return summary, [dict(assembly=source['assembly'], selection_reason=reason, windows=count, query_bp=bases)
                     for reason, (count, bases) in sorted(reasons.items())]


def one(specification):
    started = time.monotonic()
    spec, out, settings = specification
    out = Path(out)
    fasta, coverage = Path(spec['fasta']), Path(spec['coverage'])
    identity = dict(spec, fasta_sha256=digest(fasta), coverage_sha256=digest(coverage), settings=settings,
                    coverage_method='retained production samtools coverage defaults; no remapping or filtering')
    directory = out / 'assemblies' / spec['assembly']
    if directory.exists():
        previous = json.loads((directory / 'preflight_manifest.json').read_text())
        if previous['inputs'] != identity:
            raise ValueError(f"Preflight inputs/settings changed for {spec['assembly']}; choose a new --out directory")
        for relative, checksum in previous['outputs'].items():
            if digest(directory / relative) != checksum:
                raise ValueError('Changed retained preflight output: ' + str(directory / relative))
        reused = True
    else:
        with tempfile.TemporaryDirectory(prefix='.preflight-', dir=out / 'assemblies') as temporary:
            directory_tmp = Path(temporary) / 'result'
            directory_tmp.mkdir()
            metrics = directory_tmp / 'metrics'
            metrics.mkdir()
            prepare(SimpleNamespace(fasta=fasta, assembly=spec['assembly'], stage='post_fcs'), metrics)
            lengths = {row['sequence']: int(row['length_bp']) for row in table(metrics / 'sequence_metrics.tsv')}
            # Validate the complete reference dictionary; never accept a scaffold or subset table.
            subset_coverage(coverage, lengths, set(lengths), directory_tmp / 'coverage.tsv')
            sample(SimpleNamespace(fasta=fasta, metrics=metrics, coverage=directory_tmp / 'coverage.tsv',
                                   out=directory_tmp / 'queries', **settings))
            if digest(coverage) != identity['coverage_sha256']:
                raise ValueError('Coverage table changed during preparation: ' + str(coverage))
            if json.loads((metrics / 'source.json').read_text())['source_sha256'] != identity['fasta_sha256']:
                raise ValueError('Source FASTA changed before preparation: ' + str(fasta))
            files = sorted(path for path in directory_tmp.rglob('*') if path.is_file())
            manifest = dict(inputs=identity, outputs={str(path.relative_to(directory_tmp).as_posix()): digest(path)
                                                     for path in files})
            (directory_tmp / 'preflight_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
            directory_tmp.rename(directory)
        reused = False
    summary, reasons = census(directory)
    summary.update(preparation_seconds=round(time.monotonic() - started, 3), reused=reused)
    print(f"{spec['assembly']}: {summary['protein_query_bp']:,} bp "
          f"({summary['protein_query_pct']:.2f}%), {summary['protein_chunks']} protein chunks; reused={reused}", flush=True)
    return summary, reasons


def run(args):
    project = args.project.resolve(strict=True)
    trace = (args.trace or project / 'genome_assembly/pipeline/pipeline_trace.txt').resolve(strict=True)
    specs = discover(project, trace, args.expected_assemblies)
    output = args.out.resolve()
    (output / 'assemblies').mkdir(parents=True, exist_ok=True)
    settings = dict(window_bp=args.window_bp, spacing_bp=args.spacing_bp, chunk_bp=5000000,
                    short_contig_bp=args.short_contig_bp, gc_delta=args.gc_delta, depth_fold=args.depth_fold)
    # Do not leave a prior PASS looking current if a resubmission fails midway.
    (output / 'preflight.json').write_text(json.dumps(dict(status='RUNNING', settings=settings,
                                                         alignment_jobs_launched=0), indent=2) + '\n')
    started = time.monotonic()
    work = [(spec, str(output), settings) for spec in specs]
    if args.workers == 1:
        completed = [one(spec) for spec in work]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            completed = list(executor.map(one, work))
    summary = [row for row, _ in completed]
    reasons = [row for _, rows in completed for row in rows]
    write_table(output / 'summary.tsv', summary, list(summary[0]))
    write_table(output / 'selection_reasons.tsv', reasons,
                ['assembly', 'selection_reason', 'windows', 'query_bp'])
    write_table(output / 'inputs.tsv', specs, ['assembly', 'fasta', 'coverage'])
    totals = dict(assemblies=len(summary), source_bp=sum(row['source_bp'] for row in summary),
                  protein_query_bp=sum(row['protein_query_bp'] for row in summary),
                  protein_chunks=sum(row['protein_chunks'] for row in summary),
                  protein_windows=sum(row['protein_windows'] for row in summary))
    totals['protein_query_pct'] = 100 * totals['protein_query_bp'] / totals['source_bp']
    metadata = dict(status='PASS', project=str(project), trace=str(trace), settings=settings, totals=totals,
                    seconds=round(time.monotonic() - started, 3),
                    alignment_jobs_launched=0, nucleotide_query_volume='unknown until protein hits are available',
                    coverage_method='retained production samtools coverage defaults; may differ from filtered BAM coverage')
    (output / 'preflight.json').write_text(json.dumps(metadata, indent=2) + '\n')
    lines = ['# Whole-cohort staged preflight', '',
             '**Preparation only. No searches, database downloads, read mappings or removals were run.**', '',
             '| Assembly | Source Mb | Protein-query Mb | Searched % | Protein chunks |',
             '| --- | ---: | ---: | ---: | ---: |']
    lines.extend(f"| {row['assembly']} | {row['source_bp']/1e6:.3f} | {row['protein_query_bp']/1e6:.3f} | "
                 f"{row['protein_query_pct']:.2f} | {row['protein_chunks']} |" for row in summary)
    lines += ['', f"Total: {totals['protein_query_bp']:,} query bp in {totals['protein_chunks']} chunks "
              f"({totals['protein_query_pct']:.2f}% of {totals['source_bp']:,} source bp).", '',
              '[Per-assembly summary](summary.tsv) | [Selection reasons](selection_reasons.tsv) | [Source paths](inputs.tsv)', '',
              'Selection reasons are mutually exclusive combinations: summing their query_bp reproduces the total. '
              'Nucleotide volume is unknown until protein hits can be classified. No runtime or memory guarantee follows from query bp alone.', '',
              'Coverage uses retained production tables, as in the subset pilot. Filtered BAM coverage in a later pipeline run '
              'can change the depth baseline and selected windows. Per-assembly manifests bind prepared queries to exact FASTA/coverage '
              'checksums and settings. Resubmission verifies those outputs before reuse.', '']
    (output / 'README.md').write_text('\n'.join(lines))
    archive = output.with_name(output.name + '.review.tar.gz')
    with tarfile.open(archive, 'w:gz') as handle:
        for filename in ['README.md', 'summary.tsv', 'selection_reasons.tsv', 'inputs.tsv', 'preflight.json']:
            handle.add(output / filename, arcname=filename)
        for spec in specs:
            for filename in ['queries/source.json', 'metrics/source.json', 'preflight_manifest.json']:
                relative = f"assemblies/{spec['assembly']}/{filename}"
                handle.add(output / relative, arcname=relative)
    print(f"TOTAL: {totals['protein_chunks']} protein chunks; no searches launched.\nReview archive: {archive}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, default=Path.cwd())
    parser.add_argument('--trace', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--expected-assemblies', type=int, default=10)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--window-bp', type=int, default=10000)
    parser.add_argument('--spacing-bp', type=int, default=1000000)
    parser.add_argument('--short-contig-bp', type=int, default=100000)
    parser.add_argument('--gc-delta', type=float, default=.15)
    parser.add_argument('--depth-fold', type=float, default=4)
    args = parser.parse_args()
    if args.workers < 1 or args.expected_assemblies < 1:
        parser.error('workers and expected-assemblies must be positive')
    run(args)


if __name__ == '__main__':
    main()
