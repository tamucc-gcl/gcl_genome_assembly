"""Attach final assembly scope to regional naming groups without reassigning sequence."""
import argparse
import csv
import json
import re
from pathlib import Path


def chromosomes(group):
    name = group.removeprefix('composite:')
    parts = name.split('+')
    if all(re.fullmatch(r'chr\d+(?:_\d+)?', p) for p in parts):
        return [p.split('_')[0] for p in parts]
    return []


def prepare(sharing, eligibility, taxid, out):
    with open(sharing, encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assemblies = {r['id']: r for r in eligibility['assemblies'] if str(r['taxid']) == str(taxid)}
    haplotypes = {r['haplotype'] for r in rows if r['scope'] == 'haplotype'}
    if not haplotypes or haplotypes - assemblies.keys():
        raise ValueError('Regional haplotypes and species eligibility disagree')
    reported = {(r['haplotype'], r['chromosome']) for r in rows if r['scope'] == 'chromosome_haplotype'}
    expected_composites = {(hap, 'composite:'+name) for hap in haplotypes
        if assemblies[hap]['chromosome_scale'] for name in assemblies[hap]['reference_contigs'] if '+' in name}
    groups = sorted(reported | expected_composites)
    context = []
    for hap, group in groups:
        assembly = assemblies[hap]
        composite = group.startswith('composite:')
        large = bool(assembly['chromosome_scale']) and group.removeprefix('composite:') in assembly['reference_contigs']
        context.append(dict(haplotype=hap, chromosome=group,
            assembly_chromosome_scale='yes' if assembly['chromosome_scale'] else 'no',
            chromosome_scale_composite='yes' if composite and large else 'no',
            graph_group_present='yes' if (hap, group) in reported else 'no',
            kind='composite' if composite else 'named' if chromosomes(group) else 'unplaced'))
    chroms = sorted({c for r in context for c in chromosomes(r['chromosome'])}, key=lambda c: int(c[3:]))
    presence = []
    for hap in sorted(haplotypes):
        own = [r for r in context if r['haplotype'] == hap]
        for chrom in chroms:
            separate = any(r['chromosome'] == chrom for r in own)
            composites = [r for r in own if r['kind'] == 'composite' and chrom in chromosomes(r['chromosome'])]
            large = [r['chromosome'] for r in composites if r['chromosome_scale_composite'] == 'yes']
            small = [r['chromosome'] for r in composites if r['chromosome_scale_composite'] == 'no']
            missing = [r['chromosome'] for r in composites if r['graph_group_present'] == 'no']
            status = ('separate_and_composite' if separate else 'composite_only') if large else (
                'separate_group' if separate else 'other_composite_only' if small else 'no_named_group')
            if missing:
                status = 'separate_and_composite_missing' if separate else 'composite_missing_from_graph'
            presence.append(dict(haplotype=hap, chromosome=chrom, separate_group='yes' if separate else 'no',
                status=status, chromosome_scale_composites=';'.join(large), other_composites=';'.join(small),
                graph_missing_composites=';'.join(missing)))
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for name, records in [('regional_groups.tsv', context), ('chromosome_presence.tsv', presence)]:
        with (out/name).open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(records)
    return context, presence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('sharing', 'eligibility', 'taxid', 'out'):
        parser.add_argument('--'+key, required=True)
    args = parser.parse_args()
    prepare(args.sharing, json.loads(Path(args.eligibility).read_text(encoding='utf-8')), args.taxid, args.out)


if __name__ == '__main__':
    main()
