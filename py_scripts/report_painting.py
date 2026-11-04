"""Project final mapping-chain spans onto query scaffold coordinates for display.

This is reporting evidence, not a chimera detector or a cut recommendation.
Conflicting chromosome labels are explicitly ambiguous; uncovered bases unassigned.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import re


def dictionary(path):
    result = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip().split('\t')
            if fields[0] in result:
                raise ValueError('Duplicate FAI sequence: '+fields[0])
            result[fields[0]] = int(fields[1])
    return result


def chromosome(name):
    # Final harmonization names use chr<label>_<piece>. Composite reference
    # sequences do not provide an internal identity and cannot paint a query.
    if '+' in name or not name.startswith('chr'):
        return None
    return re.sub(r'_\d+$', '', name)


def partition(length, spans):
    events = defaultdict(Counter)
    events[0];events[length]
    for lo, hi, label in spans:
        if not 0 <= lo < hi <= length:
            raise ValueError('Painting interval outside source sequence')
        events[lo][label] += 1;events[hi][label] -= 1
    active = Counter();segments = []
    previous = 0
    for position in sorted(events):
        labels = [key for key, count in active.items() if count > 0]
        label = labels[0] if len(labels) == 1 else 'Ambiguous' if labels else 'Unassigned'
        if position > previous:
            if segments and segments[-1][2] == label:
                segments[-1] = (segments[-1][0], position, label)
            else:
                segments.append((previous, position, label))
        active.update(events[position])
        previous = position
    return segments


def paint(query_fai, ref_fai, paf, output, self_reference=False, min_mapq=5, min_length=150000):
    query = dictionary(query_fai);reference = dictionary(ref_fai)
    spans = defaultdict(list);counts = Counter()
    if self_reference:
        if query != reference:
            raise ValueError('Self-reference dictionaries differ')
        for name, length in query.items():
            label = chromosome(name)
            if label:
                spans[name].append((0, length, label))
    else:
        opener = gzip.open if str(paf).endswith('.gz') else open
        with opener(paf, 'rt') as handle:
            for line in handle:
                if not line.strip():continue
                f = line.rstrip().split('\t');counts['chains_seen'] += 1
                if len(f) < 12:raise ValueError('Malformed painting PAF')
                q, qlen, lo, hi, strand, target, tlen, tlo, thi, matches, aln, mapq = f[:12]
                qlen, lo, hi, tlen, tlo, thi, matches, aln, mapq = map(int, [qlen,lo,hi,tlen,tlo,thi,matches,aln,mapq])
                if q not in query or target not in reference or query[q] != qlen or reference[target] != tlen:
                    raise ValueError('PAF does not match final sequence dictionaries')
                if not (0 <= lo < hi <= qlen and 0 <= tlo < thi <= tlen and strand in ['+','-']):
                    raise ValueError('Invalid PAF coordinates')
                if mapq == 255 or mapq < min_mapq or aln < min_length or 'tp:A:S' in f[12:]:
                    counts['chains_filtered'] += 1;continue
                label = chromosome(target)
                if label is None:
                    counts['unresolved_reference_label'] += 1;continue
                spans[q].append((lo, hi, label));counts['chains_used'] += 1
    with open(output, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sequence','length','start_0based','end_exclusive','chromosome','basis'])
        for name, length in query.items():
            for lo, hi, label in partition(length, spans[name]):
                writer.writerow([name,length,lo,hi,label,'reference_identity' if self_reference else 'primary_mapping_chain'])
    measurements={key:None if self_reference else counts[key] for key in
                  ['chains_seen','chains_filtered','chains_used','unresolved_reference_label']}
    return dict(measurements, sequences=len(query), bases=sum(query.values()), self_reference=self_reference,
                min_mapq=min_mapq,min_chain_bp=min_length,coordinate_stage='final',
                interpretation='Mapping-chain spans, not base-level coverage or validation of joins')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['query-fai','ref-fai','paf','output','metadata']:p.add_argument('--'+key,required=True)
    p.add_argument('--self-reference',action='store_true')
    p.add_argument('--reference-id');p.add_argument('--query-id')
    p.add_argument('--min-mapq',type=int,default=5);p.add_argument('--min-length',type=int,default=150000)
    a=p.parse_args()
    result=paint(a.query_fai,a.ref_fai,a.paf,a.output,a.self_reference,a.min_mapq,a.min_length)
    result.update(reference_id=a.reference_id,query_id=a.query_id)
    Path(a.metadata).write_text(json.dumps(result,indent=2)+'\n')
