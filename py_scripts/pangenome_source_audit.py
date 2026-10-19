"""Trace source sequences through cached Cactus mappings and FULL/CLIP paths.

Uses the actual input FAI dictionaries, not chromosome names guessed from plots.
Mapping coverage is measured only for sequences absent under their source name
from both graph stages. It does not certify that their bases are absent elsewhere.
"""
import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path


def read_table(path):
    with open(path, encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_table(path, rows, columns):
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def union_bp(intervals):
    total, end = 0, -1
    for lo, hi in sorted(intervals):
        total += max(0, hi-max(lo, end))
        end = max(end, hi)
    return total


def query_identity(name, aliases):
    # Cactus v3.2.1 preprocessing uses id=<seqfile sample>|<FASTA ID>.
    if name.startswith('id=') and '|' in name:
        sample, seq = name[3:].split('|', 1)
        return aliases.get(sample), seq
    # PanSN names are used by exported graph paths (sample#hap#sequence).
    tokens = name.split('#', 2)
    if len(tokens) == 3:
        return aliases.get('#'.join(tokens[:2])), tokens[2]
    return None, name


def audit(specs, identities, paths, evidence, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    ledger = {r['id']: r for r in read_table(identities)}
    if not specs or {r['id'] for r in specs} != ledger.keys() or len(specs) != len(ledger):
        raise ValueError('Source audit and graph biological identities disagree')
    aliases, inputs, large, refs = {}, {}, set(), set()
    for spec in specs:
        aid = spec['id']
        row = ledger[aid]
        name = row['graph_name']
        base, dot, hap = name.rpartition('.')
        pansn = base+'#'+hap if dot and hap.isdigit() else name+'#0'
        for alias in (name, pansn):
            if alias in aliases and aliases[alias] != aid:
                raise ValueError('Ambiguous source identity alias')
            aliases[alias] = aid
        if row['is_reference'].lower() == 'true':
            refs.add(aid)
        lengths = {}
        with open(spec['fai'], encoding='utf-8') as handle:
            for line in handle:
                fields = line.rstrip().split('\t')
                if len(fields) < 2 or fields[0] in lengths or int(fields[1]) <= 0:
                    raise ValueError('Invalid or duplicate source FAI sequence')
                lengths[fields[0]] = int(fields[1])
        if len(lengths) != int(row['sequence_count']) or sum(lengths.values()) != int(row['bases']):
            raise ValueError('Source FAI does not reconcile with audited FASTA: '+aid)
        expected = set(spec['chromosomes'])
        if expected - lengths.keys():
            raise ValueError('Expected chromosome missing from source FAI: '+aid)
        large.update((aid, seq) for seq in expected)
        inputs.update({(aid, seq): n for seq, n in lengths.items()})
    graph = {}
    for row in read_table(paths):
        key = row['assembly'], row['sequence']
        if key in graph or key not in inputs:
            raise ValueError('Unknown or duplicate graph source path: '+str(key))
        full = int(row['full_W_coordinate_span_bp']) if row['full_W_coordinate_span_bp'] else None
        clip = int(row['clip_retained_path_bp'])
        if full is not None and not 0 <= full <= inputs[key]:
            raise ValueError('FULL source span exceeds input length')
        if clip < 0:
            raise ValueError('Negative CLIP retained path length')
        graph[key] = full, clip
    missing = {key for key in inputs if key not in graph or not any(v for v in graph[key])}
    measurements, inventory, unknown, excerpts = {}, [], Counter(), []
    seen_kinds = set()
    for item in evidence:
        kind, path = item['kind'], Path(item['path'])
        if kind not in ('gaf', 'paf') or kind in seen_kinds:
            raise ValueError('Duplicate or unsupported Cactus mapping evidence')
        seen_kinds.add(kind)
        counts, intervals, highq = Counter(), defaultdict(list), Counter()
        examples = Counter()
        opener = gzip.open if path.suffix == '.gz' else open
        total, synthetic = 0, 0
        with opener(path, 'rt', encoding='utf-8') as handle:
            for line in handle:
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip().split('\t')
                if len(fields) < 12:
                    raise ValueError('Malformed '+kind+' record')
                total += 1
                if fields[0].startswith(('id=_MINIGRAPH_|','_MINIGRAPH_#')):
                    synthetic += 1
                    continue
                aid, seq = query_identity(fields[0], aliases)
                if not aid:
                    unknown[fields[0].split('|', 1)[0].split('#', 2)[0]] += 1
                    continue
                key = aid, seq
                if key not in inputs:
                    raise ValueError('Mapping query absent from source FAI: '+str(key))
                qlen, lo, hi, mapq = int(fields[1]), int(fields[2]), int(fields[3]), int(fields[11])
                if qlen != inputs[key] or not 0 <= lo <= hi <= qlen or not 0 <= mapq <= 255:
                    raise ValueError('Mapping query coordinates disagree with input: '+str(key))
                if key not in missing:
                    continue
                counts[key] += 1
                intervals[key].append((lo, hi))
                highq[key] += 20 <= mapq < 255
                if examples[key] < 20:
                    excerpts.append(dict(kind=kind, assembly=aid, sequence=seq,
                        query_start=lo, query_end=hi, mapq=mapq, target=fields[5][:2000]))
                    examples[key] += 1
        measurements[kind] = dict(counts=counts, intervals=intervals, highq=highq)
        inventory.append(dict(kind=kind, file=path.name, bytes=path.stat().st_size, records=total,
                              synthetic_query_records=synthetic))
    rows = []
    for key, length in sorted(inputs.items()):
        aid, seq = key
        full, clip = graph.get(key, (None, 0))
        if key not in missing:
            status = 'represented_in_clip' if clip else 'full_only'
        elif ('gaf' in measurements and measurements['gaf']['counts'][key]
              and 'paf' in measurements and not measurements['paf']['counts'][key]):
            status = 'gaf_present_but_no_paf_query_records'
        elif 'gaf' in measurements and measurements['gaf']['counts'][key]:
            status = 'mapped_but_absent_from_full_and_clip'
        elif 'paf' in measurements and measurements['paf']['counts'][key]:
            status = 'paf_present_but_absent_from_full_and_clip'
        elif aid in refs:
            status = 'reference_path_missing_mapping_not_diagnostic'
        elif not seen_kinds:
            status = 'missing_path_mapping_unavailable'
        else:
            status = 'no_query_records_in_available_mappings'
        row = dict(assembly=aid, sequence=seq, source_bp=length,
                   chromosome_scale='yes' if key in large else 'no',
                   full_W_coordinate_span_bp=full if full is not None else '',
                   clip_retained_path_bp=clip, status=status)
        for kind in ('gaf','paf'):
            measured = measurements.get(kind)
            assess = measured is not None and key in missing
            row[kind+'_records'] = measured['counts'][key] if assess else ''
            row[kind+'_query_covered_bp'] = union_bp(measured['intervals'][key]) if assess else ''
            row[kind+'_mapq20_records'] = measured['highq'][key] if assess else ''
        rows.append(row)
    problems = [r for r in rows if r['chromosome_scale']=='yes' and r['status']!='represented_in_clip']
    result = dict(status='REVIEW' if problems or unknown else 'PASS',
        input_sequences=len(inputs), source_names_missing_from_both_graphs=len(missing),
        chromosome_scale_paths_requiring_review=len(problems),
        unknown_mapping_name_prefixes=dict(unknown), mapping_evidence=inventory)
    write_table(out/'source_sequences.tsv', rows, list(rows[0]))
    write_table(out/'missing_source_paths.tsv', [r for r in rows if (r['assembly'],r['sequence']) in missing], list(rows[0]))
    write_table(out/'mapping_examples.tsv', excerpts,
                ['kind','assembly','sequence','query_start','query_end','mapq','target'])
    (out/'source_audit.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    lines = ['# Source sequence representation', '',
        f"Status: **{result['status']}**. {len(inputs):,} input sequences checked; {len(problems)} chromosome-scale paths require review.", '',
        'All input FAI dictionaries reconcile with the graph input audit. Scope comes from the final assembly eligibility audit.', '',
        '| Assembly | Sequence | Source bp | GAF covered bp | PAF covered bp | Finding |', '|---|---|---:|---:|---:|---|']
    labels = dict(full_only='Present in FULL; no retained CLIP path',
        mapped_but_absent_from_full_and_clip='Mapped by minigraph; source name absent from FULL and CLIP',
        gaf_present_but_no_paf_query_records='GAF query present; no query in exported PAF or FULL/CLIP',
        paf_present_but_absent_from_full_and_clip='Present in mapping PAF; source name absent from FULL and CLIP',
        reference_path_missing_mapping_not_diagnostic='Reference path missing; query mapping is not diagnostic',
        missing_path_mapping_unavailable='Source name missing; mapping evidence unavailable',
        no_query_records_in_available_mappings='Source name missing; no query records in available mappings')
    def show(value):
        return f'{value:,}' if value != '' else 'Not measured'
    lines += [f"| {r['assembly']} | {r['sequence']} | {r['source_bp']:,} | {show(r['gaf_query_covered_bp'])} | {show(r['paf_query_covered_bp'])} | {labels[r['status']]} |" for r in problems]
    if not problems:
        lines += ['All expected chromosome-scale source paths have retained CLIP paths.']
    lines += ['', '## Interpretation', '',
        'A missing source name is a representation discrepancy, not proof of a biological deletion or total loss of its bases. '
        'Non-reference query records present upstream but absent from FULL localize the discrepancy after mapping; '
        'they do not distinguish chromosome partitioning from subsequent alignment or graph editing. '
        'No query record means none in the supplied files, not that no homologous sequence exists. '
        'Reference mapping can be extracted directly from graph tags rather than appear as a minigraph query.', '',
        'Coverage is the union of query intervals, never a sum of overlapping alignments. MAPQ 255 is unavailable, '
        'not high quality. GAF MAPQ counts are descriptive and do not reimplement Cactus filtering. '
        'Only missing source paths have mapping measurements; blank cells for other paths mean not measured. '
        'Missing mapping files are reported explicitly, never measured zero.', '',
        '[All source sequences](source_sequences.tsv) · [Missing source names](missing_source_paths.tsv) · '
        '[Mapping examples (first 20 per sequence/file; targets truncated to 2,000 characters)](mapping_examples.tsv) · [Audit and file inventory](source_audit.json)', '',
        'Chromosome-assignment logs from older cached builds may have been removed with construction scratch. '
        'This audit does not reconstruct their contents or claim a specific partitioning decision without them.', '']
    (out/'source_report.md').write_text('\n'.join(lines), encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('input','identities','paths','out'):
        parser.add_argument('--'+key, required=True)
    args = parser.parse_args()
    spec = json.loads(Path(args.input).read_text(encoding='utf-8'))
    audit(spec['assemblies'], args.identities, args.paths, spec['evidence'], args.out)


if __name__ == '__main__':
    main()
