#!/usr/bin/env python3
"""Join exact-screen FASTA metrics, read coverage and reported FCS classifications.

No similarity searches or contamination selection thresholds are introduced here.
Inconclusive FCS results remain inconclusive; segment disagreements remain mixed.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re

from contamination_review import records


def read_taxonomy(path):
    lines = Path(path).read_text().splitlines()
    header = next((line for line in lines if line.startswith('#seq-id\t')), None)
    if header is None:
        raise ValueError('Missing FCS taxonomy header: ' + str(path))
    return list(csv.DictReader([header.lstrip('#')] + [x for x in lines if x and not x.startswith('#')], delimiter='\t'))


def parent_sequence(identifier, lengths):
    # GX can report a segment of a segment (e.g. scaffold~a..b~~c..d).
    # Stop at the nearest exact FASTA name: coordinate-looking text can be a
    # literal part of that name. Then validate every suffix in its parent's frame.
    name=identifier
    segments=[]
    while name not in lengths:
        match = re.fullmatch(r'(.+?)(?:~~|~)(\d+)\.\.(\d+)', name)
        if not match:
            raise ValueError('FCS sequence absent from exact screened FASTA: ' + identifier)
        segments.append((int(match[2]),int(match[3])))
        name=match[1]
    parent_length=lengths[name]
    for lo,hi in reversed(segments):
        if not 1 <= lo <= hi <= parent_length:
            raise ValueError('FCS segment coordinates outside source: ' + identifier)
        parent_length=hi-lo+1
    return name


def classification(rows):
    if not rows:
        return 'unassigned'
    divisions = {x['div'] for x in rows if x.get('div') and x['div'] not in ('n/a','NA')}
    results = {x['result'] for x in rows}
    confident = all(x.startswith(('primary-div','contaminant')) for x in results)
    if len(divisions) == 1 and confident:
        return next(iter(divisions))
    if len(divisions) > 1 or (len(results) > 1 and any(x.startswith('contaminant') for x in results)):
        return 'mixed FCS classifications'
    return 'unresolved: ' + ', '.join(sorted(results))


def build(fasta, taxonomy, coverage, assembly, stage, out):
    metrics = []
    seen = set()
    for header, sequence in records(fasta):
        name = header.split()[0]
        if name in seen:
            raise ValueError('Duplicate FASTA identifier: ' + name)
        seen.add(name)
        upper = sequence.upper()
        acgt = sum(upper.count(base) for base in 'ACGT')
        metrics.append(dict(sequence=name,length_bp=len(sequence),acgt_bp=acgt,
                            gc_fraction=(upper.count('G')+upper.count('C'))/acgt if acgt else 'NA'))
    lengths = {x['sequence']:x['length_bp'] for x in metrics}
    grouped = defaultdict(list)
    tax_rows = read_taxonomy(taxonomy)
    for row in tax_rows:
        identifier = row['seq-id']
        name = parent_sequence(identifier,lengths)
        if identifier == name and int(row['seq-len']) != lengths[name]:
            raise ValueError('FCS taxonomy length differs from screened FASTA: ' + name)
        grouped[name].append(row)
    with Path(coverage).open() as handle:
        depth_rows = list(csv.DictReader(handle,delimiter='\t'))
    depths = {}
    for row in depth_rows:
        name = row['#rname']
        if name == '*':
            continue
        if name in depths or name not in lengths:
            raise ValueError('Duplicate or unknown coverage sequence: ' + name)
        if int(row['startpos']) != 1 or int(row['endpos']) != lengths[name]:
            raise ValueError('Coverage coordinates differ from screened FASTA: ' + name)
        depths[name] = row['meandepth']
    if set(depths) != set(lengths):
        raise ValueError('Coverage must include every screened contig, including zero-depth contigs')
    for row in metrics:
        source = grouped[row['sequence']]
        row.update(mean_depth=depths[row['sequence']],plot_taxon=classification(source),
                   fcs_results=';'.join(sorted({x['result'] for x in source})),
                   fcs_divisions=';'.join(sorted({x['div'] for x in source})),
                   classification_source='FCS-GX',assembly=assembly,coordinate_stage=stage)
    out = Path(out)
    with out.open('w',newline='') as handle:
        writer = csv.DictWriter(handle,list(metrics[0]) if metrics else
                                ['sequence','length_bp','gc_fraction','mean_depth','plot_taxon'],delimiter='\t')
        writer.writeheader()
        writer.writerows(metrics)
    def sha(path):
        result = hashlib.sha256()
        with Path(path).open('rb') as handle:
            for chunk in iter(lambda:handle.read(8*1024*1024),b''):
                result.update(chunk)
        return result.hexdigest()
    out.with_suffix('.json').write_text(json.dumps(dict(assembly=assembly,coordinate_stage=stage,
        classification_source='FCS-GX taxonomy report; no BLAST or sampled assignment',
        source_sha256=sha(fasta),taxonomy_sha256=sha(taxonomy),coverage_sha256=sha(coverage),
        contigs=len(metrics),taxonomy_rows=len(tax_rows)),indent=2)+'\n')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('fasta','taxonomy','coverage','out'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--assembly',required=True)
    parser.add_argument('--stage',required=True)
    args=parser.parse_args()
    build(args.fasta,args.taxonomy,args.coverage,args.assembly,args.stage,args.out)
