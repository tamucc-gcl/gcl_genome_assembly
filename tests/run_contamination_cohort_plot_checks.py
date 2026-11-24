"""Render real GC/coverage figures and verify missing/zero-depth accounting."""
import argparse
import csv
from pathlib import Path
import subprocess


def check(rscript, root, label, records, expected):
    out = root / label
    out.mkdir(parents=True, exist_ok=True)
    source = out / 'evidence.tsv'
    fields = ['cohort_assembly', 'sequence', 'length_bp', 'gc_fraction', 'mean_depth', 'plot_taxon']
    with source.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(fields)
        writer.writerows(records)
    result = subprocess.run([rscript, str(Path(__file__).resolve().parents[1] / 'r_scripts/plot_contamination_cohort.R'),
                             str(source), str(out)], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    with (out / 'cohort_gc_coverage_scope.tsv').open() as handle:
        scope = {row['assembly']: row for row in csv.DictReader(handle, delimiter='\t')}
    for assembly, counts in expected.items():
        for name, value in counts.items():
            assert int(scope[assembly][name]) == value, (label, assembly, name, scope[assembly])
    assert (out / 'cohort_gc_coverage.png').read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    assert (out / 'cohort_gc_coverage.pdf').read_bytes().startswith(b'%PDF')
    print('PASS', label)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rscript')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1] / 'comparisons/contamination-cohort-plot-validation'
    check(args.rscript, root, 'mixed', [
        ['A', 'zero', 1000, .4, 0, 'Chordata'], ['A', 'depth', 500, .5, 20, 'Chordata'],
        ['A', 'no_gc', 100, 'NA', 10, 'no-hit'], ['B', 'no_depth', 100, .6, 'NA', 'unresolved']],
        {'A': dict(contigs=3, plotted=2, zero_depth=1, omitted=1), 'B': dict(contigs=1, plotted=0, omitted=1)})
    check(args.rscript, root, 'no_depth', [['A', 'a', 100, .5, 'NA', 'no-hit']],
          {'A': dict(contigs=1, plotted=0, missing_depth=1)})
    check(args.rscript, root, 'single', [['A', 'a', 100, .5, 20, 'Chordata']],
          {'A': dict(contigs=1, plotted=1, omitted=0)})
