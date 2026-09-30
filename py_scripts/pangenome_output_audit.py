#!/usr/bin/env python3
"""Lightweight artifact checks; no graph traversal or biological validation."""
import argparse
import csv
import gzip
import json
from pathlib import Path

REQUIRED = ('biological_gbz', 'biological_gfa', 'biological_odgi', 'haplotype_index',
            'variant_gbz', 'variant_gfa', 'variants_standard', 'variants_standard_index',
            'full_qc_gfa', 'clipping_stats')


def audit(roles):
    rows, errors = [], []
    for role in dict.fromkeys((*REQUIRED, *roles)):
        name = roles.get(role, '')
        path = Path(name) if name else None
        status = 'present' if path and path.is_file() and path.stat().st_size else 'missing_or_empty'
        detail = ''
        if status == 'present' and role == 'variants_standard':
            try:
                with gzip.open(path, 'rt') as handle:
                    if not handle.readline().startswith('##fileformat=VCFv'):
                        raise ValueError('missing VCF fileformat header')
                    for line in handle:
                        if line.startswith('#CHROM\t'):
                            fields = line.rstrip('\n').split('\t')
                            if fields[:8] != ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']:
                                raise ValueError('invalid VCF column header')
                            detail = 'VCF header readable; %d sample columns' % max(0, len(fields)-9)
                            break
                        if not line.startswith('#'):
                            raise ValueError('variant before column header')
                    else:
                        raise ValueError('missing VCF column header')
            except (OSError, EOFError, UnicodeError, ValueError) as exc:
                status, detail = 'invalid', str(exc)
        graph = 'gref_clip' if role.startswith('variant') else 'full_support' if role.startswith('full_') or role == 'clipping_stats' else 'clip'
        rows.append(dict(role=role, graph=graph, file=name, required=role in REQUIRED,
                         status=status, bytes=path.stat().st_size if path and path.is_file() else 0, detail=detail))
        if role in REQUIRED and status != 'present':
            errors.append(role+': '+status+(' ('+detail+')' if detail else ''))
    return rows, errors


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--roles', required=True)
    a = p.parse_args()
    rows, errors = audit(json.loads(Path(a.roles).read_text()))
    with open('pangenome_artifact_checks.tsv', 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    if errors:
        raise SystemExit('Graph output audit failed; retain cached Cactus outputs. '+ '; '.join(errors))


if __name__ == '__main__':
    main()
