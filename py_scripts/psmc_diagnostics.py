"""Export screening diagnostics for an individual PSMC input and observed fit."""
import argparse
import csv
from itertools import accumulate
import json
from pathlib import Path

from psmc_results import parse, pattern_groups, read_fasta


def write_tsv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def input_metrics(input_dir, qc, record):
    category = 'Input quality'
    bins = qc['bins']
    usable = bins['K'] + bins['T']
    bin_size = qc['bin_size']
    with (input_dir / 'contig_qc.tsv').open(encoding='utf-8') as handle:
        lengths = {r['contig']: int(r['length_bp']) for r in csv.DictReader(handle, delimiter='\t')}
    callable_span = 0
    for contig, sequence in read_fasta(input_dir / 'input.psmcfa'):
        length = lengths[contig]
        callable_span += (sequence.count('K') + sequence.count('T')) * bin_size
        if sequence[-1] != 'N':
            callable_span -= (-length) % bin_size

    record(category, 'Fit eligibility', qc['status'], status=qc['reason'] or 'reported')
    record(category, 'Callable bases before binning', qc['callable_bp'], 'bp')
    record(category, 'Callable fraction of selected non-gap scope', qc['callable_fraction'], 'fraction',
           'low callability' if qc['callable_fraction'] < 0.5 else 'reported')
    record(category, 'Heterozygosity per callable base', qc['heterozygosity'], 'fraction',
           'review heterozygosity' if qc['heterozygosity'] is not None
           and not 0.0001 <= qc['heterozygosity'] <= 0.02 else 'reported')
    record(category, 'Retained input sequences', qc['segments'], 'sequences')
    record(category, 'Retained callable bases before binning', qc['psmc_retained_callable_bp'], 'bp')
    record(category, 'Usable bins', usable, 'bins')
    record(category, 'Heterozygous bins', bins['K'], 'bins')
    record(category, 'Masked fraction of retained bins',
           bins['N'] / sum(bins.values()) if sum(bins.values()) else None, 'fraction',
           'high masking' if bins['N'] and bins['N'] >= usable else 'reported')
    record(category, 'Fraction of usable bins containing heterozygotes',
           bins['K'] / usable if usable else None, 'fraction')
    record(category, 'Reference fraction spanned by usable bins',
           callable_span / sum(lengths.values()) if sum(lengths.values()) else None, 'fraction')
    if qc['mode'] == 'assemblies':
        outside = qc.get('mask_audit_totals', {}).get('outside_confident_callable_bp')
        record(category, 'Callable bases outside Dipcall confident BED', outside, 'bp',
               'review mask leakage' if outside else 'pass')


def observed_metrics(observed, qc, out, iteration, mutation_rate, generation_time, record):
    fit = parse(observed)
    blocks = sorted(fit['blocks'], key=lambda block: block['iteration'])
    chosen = next(block for block in blocks if block['iteration'] == iteration)
    groups = pattern_groups(fit['pattern'])
    bounds = [0, *accumulate(groups)]
    theta_divisor = qc['bin_size'] * fit['skip']
    scaled = mutation_rate is not None
    time_unit = 'years' if scaled else 'scaled divergence'
    population_unit = 'Ne' if scaled else 'scaled size'

    theta = chosen['theta'] / theta_divisor
    time_factor = theta * generation_time / (2 * mutation_rate) if scaled else theta
    size_factor = theta / (4 * mutation_rate) if scaled else theta
    intervals = []
    for group, (start, stop) in enumerate(zip(bounds, bounds[1:]), 1):
        rows = chosen['rows']
        count = sum(row[3] for row in rows[start:stop])
        shift = chosen.get('dt', 0)
        intervals.append(dict(
            parameter_group=group, iteration=iteration, first_atomic_interval=start,
            last_atomic_interval=stop - 1, time_start=(rows[start][1] + shift) * time_factor,
            time_end=(rows[stop][1] + shift) * time_factor if stop < len(rows) else None,
            time_unit=time_unit, relative_size=rows[start][2], population=rows[start][2] * size_factor,
            population_unit=population_unit, expected_recombinations=count,
            support='insufficient' if count < 10 else 'caution' if count < 20 else 'reasonable',
            recent_peak='no'))

    # Inspect the first five groups, excluding endpoints without two neighbors.
    peaks = []
    for younger, middle, older in list(zip(intervals, intervals[1:], intervals[2:]))[:4]:
        if middle['population'] >= 2 * max(younger['population'], older['population']):
            middle['recent_peak'] = 'yes'
            peaks.append(middle['parameter_group'])
    write_tsv(out / 'psmc_intervals.tsv', intervals)

    support = 'Recombination support'
    counts = [interval['expected_recombinations'] for interval in intervals]
    record(support, 'Diagnostic iteration', iteration)
    record(support, 'Minimum expected recombinations', min(counts), 'events', intervals[counts.index(min(counts))]['support'])
    record(support, 'Intervals below 10 events', sum(n < 10 for n in counts), 'intervals',
           'insufficient' if min(counts) < 10 else 'reported')
    record(support, 'Intervals below 20 events', sum(n < 20 for n in counts), 'intervals',
           'review' if min(counts) < 20 else 'reported')
    supported = [interval for interval in intervals if interval['support'] == 'reasonable']
    for label, interval in [('Youngest', supported[0] if supported else None),
                            ('Oldest', supported[-1] if supported else None)]:
        record('Temporal support', f'{label} interval with at least 20 events',
               interval['parameter_group'] if interval else None, 'group')
        if interval:
            record('Temporal support', f'{label} supported interval start', interval['time_start'], time_unit)
            record('Temporal support', f'{label} supported interval end', interval['time_end'], time_unit,
                   'open-ended' if interval['time_end'] is None else 'reported')
    record('Peak screening', 'Recent peak groups', ', '.join(map(str, peaks)) if peaks else
           'None' if len(intervals) >= 3 else None, status='review' if peaks else 'reported')

    history = [dict(iteration=block['iteration'], log_likelihood=block.get('log_likelihood'),
                    theta_per_base=block['theta'] / theta_divisor) for block in blocks]
    write_tsv(out / 'psmc_iterations.tsv', history)
    final = blocks[-1]
    # Compare the final fit to the fit five updates earlier when available.
    start = next((b for b in reversed(blocks[:-1]) if b['iteration'] <= final['iteration'] - 5), blocks[0])
    category = 'Convergence'
    record(category, 'PSMC version', fit['version'])
    record(category, 'Interval pattern', fit['pattern'])
    record(category, 'Free interval parameters', len(groups))
    record(category, 'Iterations completed', final['iteration'])
    record(category, 'Final log-likelihood', final.get('log_likelihood'))
    likelihoods = [b for b in blocks if b['iteration'] >= start['iteration'] and 'log_likelihood' in b]
    if len(likelihoods) >= 2:
        first, last = likelihoods[0], likelihoods[-1]
        record(category, 'Log-likelihood comparison rounds', f"{first['iteration']} to {last['iteration']}")
        record(category, 'Change in log-likelihood', last['log_likelihood'] - first['log_likelihood'])
    else:
        record(category, 'Change in log-likelihood', None)
    record(category, 'Size comparison rounds', f"{start['iteration']} to {final['iteration']}")
    relative_changes = [final['rows'][i][2] / start['rows'][i][2] for i in bounds[:-1]]
    theta_ratio = final['theta'] / start['theta']
    record(category, 'Maximum change in relative size', 100 * max(abs(ratio - 1) for ratio in relative_changes), '%')
    record(category, f'Maximum change in {population_unit}',
           100 * max(abs(ratio * theta_ratio - 1) for ratio in relative_changes), '%')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--observed', type=Path)
    parser.add_argument('--iteration', type=int, required=True)
    parser.add_argument('--mutation-rate', type=float)
    parser.add_argument('--generation-time', type=float)
    args = parser.parse_args()

    qc = json.loads((args.input_dir / 'input_qc.json').read_text())
    if (qc['status'] == 'ready') != (args.observed is not None):
        parser.error('An observed fit is required exactly when input is ready')
    args.out.mkdir(parents=True, exist_ok=True)
    metrics = []

    def record(category, metric, value, unit='', status='reported'):
        if value is None and status != 'open-ended':
            status = 'unavailable'
        metrics.append(dict(category=category, metric=metric, value=value, unit=unit, status=status))

    if args.observed is not None:
        observed_metrics(args.observed, qc, args.out, args.iteration, args.mutation_rate,
                         args.generation_time, record)
    input_metrics(args.input_dir, qc, record)
    write_tsv(args.out / 'psmc_qc.tsv', metrics)


if __name__ == '__main__':
    main()
