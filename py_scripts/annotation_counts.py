#!/usr/bin/env python3
"""Read-only GenBank locus accounting for reports; no annotation is modified."""
import argparse
import csv
from pathlib import Path
import Bio
from Bio import SeqIO


def audit(path, sample):
    with Path(path).open(encoding='utf-8') as handle:
        records = list(SeqIO.parse(handle, 'genbank'))
    if not records:
        raise ValueError('No GenBank records: '+str(path))
    if len({r.id for r in records}) != len(records):
        raise ValueError('Duplicate GenBank record identity')
    loci = {}
    feature_counts = dict(CDS=set(), tRNA=set(), rRNA=set())
    for record in records:
        for feature in record.features:
            if feature.type not in ('gene', 'CDS', 'tRNA', 'rRNA'):
                continue
            if feature.location is None:
                raise ValueError('Annotated feature has no location')
            parts = tuple(sorted((int(p.start), int(p.end), p.strand or 0) for p in feature.location.parts))
            if any(lo < 0 or hi > len(record.seq) or hi <= lo for lo, hi, strand in parts):
                raise ValueError('Annotated feature lies outside its record')
            key = (record.id, parts)
            locus = loci.setdefault(key, dict(record=record.id, location=str(feature.location), types=set(), genes=set(), products=set(), feature_rows=0))
            locus['types'].add(feature.type)
            locus['genes'].update(feature.qualifiers.get('gene', []))
            locus['products'].update(feature.qualifiers.get('product', []))
            locus['feature_rows'] += 1
            if feature.type in feature_counts:
                feature_counts[feature.type].add(key)
    rows = [dict(sample_id=sample, record=l['record'], location=l['location'],
                 feature_types=';'.join(sorted(l['types'])), genes=';'.join(sorted(l['genes'])),
                 products=';'.join(sorted(l['products'])), annotation_feature_rows=l['feature_rows'])
            for key, l in sorted(loci.items())]
    summary = dict(sample_id=sample, annotation_length_bp=sum(len(r.seq) for r in records),
                   annotated_loci=len(loci), protein_coding_loci=len(feature_counts['CDS']),
                   trna_loci=len(feature_counts['tRNA']), rrna_loci=len(feature_counts['rRNA']),
                   other_gene_loci=sum(not (l['types'] & {'CDS','tRNA','rRNA'}) for l in loci.values()),
                   biopython_version=Bio.__version__)
    return summary, rows


def write(path, rows):
    with Path(path).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genbank', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()
    summary, rows = audit(args.genbank, args.sample)
    write(args.prefix+'.counts.tsv', [summary])
    if rows:
        write(args.prefix+'.loci.tsv', rows)
    else:
        Path(args.prefix+'.loci.tsv').write_text('sample_id\trecord\tlocation\tfeature_types\tgenes\tproducts\tannotation_feature_rows\n')


if __name__ == '__main__':
    main()
