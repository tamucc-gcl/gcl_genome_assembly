"""Summarize Cactus clipping evidence and verification-image provenance after construction.

Cactus v3.2.1 graphmap_join exports sample-stats.tsv (graph partition, PanSN
sample/haplotype, retained path bp) and clip-stats.bed in its stats.tgz bundle.
No graph reconstruction or source-interval attribution is inferred here.
"""
import argparse
import csv
import gzip
import importlib.metadata
import json
from collections import defaultdict
from pathlib import Path, PurePosixPath
import tarfile


def write_table(path, rows, columns):
    with Path(path).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def union_bp(intervals):
    total = 0
    end = -1
    for lo, hi in sorted(intervals):
        total += max(0, hi-max(lo, end))
        end = max(end, hi)
    return total


def summarize(bundle, full_gfa, clip_stats, identities, out, image, build_versions):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with open(identities, encoding='utf-8') as handle:
        ledger = list(csv.DictReader(handle, delimiter='\t'))
    if not ledger or len({r['id'] for r in ledger}) != len(ledger):
        raise ValueError('Empty or duplicate biological identities')
    aliases = {}
    for row in ledger:
        name = row['graph_name']
        base, dot, hap = name.rpartition('.')
        pansn = base+'#'+hap if dot and hap.isdigit() else name+'#0'
        for key in (name, pansn):
            if key in aliases and aliases[key] != row['id']:
                raise ValueError('Ambiguous graph identity alias')
            aliases[key] = row['id']

    def identity(name):
        # Cactus path/sample stats use PanSN; clip BED may carry full paths.
        if name in aliases:
            return aliases[name]
        if '#' in name:
            return aliases.get('#'.join(name.split('#')[:2]))
        return None

    retained = defaultdict(int)
    by_partition = defaultdict(int)
    clip_paths = defaultdict(int)
    regions = defaultdict(list)
    unknown = defaultdict(lambda: dict(rows=0, bp=0))
    inventory = []
    required = {'sample-stats.tsv', 'path-stats.tsv', 'clip-stats.bed'}
    with tarfile.open(bundle, 'r:gz') as archive:
        members = {}
        for member in archive:
            if not member.isfile():
                continue
            key = PurePosixPath(member.name).name
            if key in members:
                raise ValueError('Duplicate clipping bundle member: '+key)
            members[key] = member
        if required - members.keys():
            raise ValueError('Cactus clipping bundle lacks: '+', '.join(sorted(required-members.keys())))
        for name, member in members.items():
            inventory.append(dict(member=name, bytes=member.size))
        for name in sorted(required):
            with archive.extractfile(members[name]) as handle:
                for raw in handle:
                    line = raw.decode('utf-8').strip()
                    if not line or line.startswith('#'):
                        continue
                    fields = line.split('\t')
                    if name in ('sample-stats.tsv', 'path-stats.tsv'):
                        if len(fields) != 3:
                            raise ValueError('Unexpected Cactus sample-stats schema')
                        partition, sample, value = fields
                        bp = int(value)
                        if bp < 0:
                            raise ValueError('Negative retained path length')
                        aid = identity(sample)
                        if aid:
                            if name == 'sample-stats.tsv':
                                retained[aid] += bp
                                by_partition[partition, aid] += bp
                            else:
                                tokens = sample.split('#', 2)
                                if len(tokens) != 3:
                                    raise ValueError('Non-PanSN Cactus path-stats name: '+sample)
                                clip_paths[aid, tokens[2].split('[',1)[0]] += bp
                        else:
                            unknown[name+':'+sample]['rows'] += 1
                            unknown[name+':'+sample]['bp'] += bp
                    else:
                        if len(fields) < 3:
                            raise ValueError('Unexpected Cactus clip BED schema')
                        path, lo, hi = fields[0], int(fields[1]), int(fields[2])
                        if lo < 0 or hi < lo:
                            raise ValueError('Invalid clipped interval')
                        aid = identity(path)
                        if aid:
                            regions[aid, path].append((lo, hi))
                        else:
                            unknown['clipped:'+path]['rows'] += 1
                            unknown['clipped:'+path]['bp'] += hi-lo
    clipped = defaultdict(int)
    for (aid, path), intervals in regions.items():
        clipped[aid] += union_bp(intervals)
    if {r['id'] for r in ledger} - retained.keys():
        raise ValueError('Cactus sample-stats lacks audited biological assemblies')
    for aid, total in retained.items():
        if sum(v for (a, s), v in clip_paths.items() if a == aid) != total:
            raise ValueError('Cactus path/sample retained counts disagree: '+aid)
    # Unresolved tool names are explicit missing evidence, never measured zero.
    records = []
    for row in ledger:
        aid, source = row['id'], int(row['bases'])
        records.append(dict(assembly=aid, source_fasta_bp=source,
            retained_path_bp=retained.get(aid, 0),
            retained_path_percent=100*retained.get(aid, 0)/source,
            reported_clip_interval_bp=clipped.get(aid, 0),
            clip_identity_status='incomplete_name_mapping' if unknown else 'mapped',
            input_minus_retained_minus_reported_clip_bp=source-retained.get(aid, 0)-clipped.get(aid, 0)))
    write_table(out/'assembly_clipping.tsv', records, list(records[0]))
    write_table(out/'partition_retained_paths.tsv',
        [dict(partition=p, assembly=a, retained_path_bp=n) for (p, a), n in sorted(by_partition.items())],
        ['partition', 'assembly', 'retained_path_bp'])
    write_table(out/'unmapped_tool_names.tsv', [dict(name=k, **v) for k, v in sorted(unknown.items())], ['name','rows','bp'])
    write_table(out/'clipping_bundle_inventory.tsv', inventory, ['member', 'bytes'])
    full = dict(nodes=0, sequence_bp=0, walk_records=0)
    full_intervals = defaultdict(list)
    with gzip.open(full_gfa, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('S\t'):
                fields = line.rstrip().split('\t')
                full['nodes'] += 1
                if fields[2] == '*':
                    length = next((int(t[5:]) for t in fields[3:] if t.startswith('LN:i:')), None)
                    if length is None:
                        raise ValueError('FULL segment lacks sequence and LN tag')
                else:
                    length = len(fields[2])
                full['sequence_bp'] += length
            elif line.startswith(('W\t', 'P\t')):
                full['walk_records'] += 1
                if line.startswith('W\t'):
                    fields = line.rstrip().split('\t',6)
                    aid = identity('#'.join(fields[1:3]))
                    if aid and fields[4].isdigit() and fields[5].isdigit():
                        lo, hi = int(fields[4]), int(fields[5])
                        if hi < lo:
                            raise ValueError('Invalid FULL walk coordinates')
                        full_intervals[aid, fields[3]].append((lo,hi))
    clip = {}
    for line in Path(clip_stats).read_text(encoding='utf-8').splitlines():
        fields = line.split('\t')
        if len(fields) == 2 and fields[0] in ('nodes','edges','length'):
            clip[fields[0]] = int(fields[1])
    if not full['nodes'] or not full['sequence_bp'] or set(clip) != {'nodes','edges','length'}:
        raise ValueError('Empty FULL graph or incomplete CLIP statistics')
    write_table(out/'full_clip_graph.tsv', [
        dict(graph='FULL_support', nodes=full['nodes'], sequence_bp=full['sequence_bp']),
        dict(graph='CLIP_biological', nodes=clip['nodes'], sequence_bp=clip['length'])],
        ['graph','nodes','sequence_bp'])
    path_rows = [dict(assembly=aid, sequence=seq,
        full_W_coordinate_span_bp=union_bp(full_intervals[aid,seq]) if (aid,seq) in full_intervals else '',
        clip_retained_path_bp=clip_paths.get((aid,seq),0))
        for aid,seq in sorted(full_intervals.keys() | clip_paths.keys())]
    write_table(out/'full_clip_sequence_paths.tsv',path_rows,
                ['assembly','sequence','full_W_coordinate_span_bp','clip_retained_path_bp'])
    lost = sorted([r for r in path_rows if r['full_W_coordinate_span_bp'] and not r['clip_retained_path_bp']],
                  key=lambda r: -r['full_W_coordinate_span_bp'])
    versions = dict(verification_image=image,
        original_build_version_record=Path(build_versions).read_text(encoding='utf-8').strip())
    for tool in ('cactus','toil'):
        try:
            versions[tool+'_package'] = importlib.metadata.version(tool)
        except importlib.metadata.PackageNotFoundError:
            versions[tool+'_package'] = 'not_available'
    try:
        from cactus.shared.version import cactus_commit
        versions['cactus_commit'] = cactus_commit
    except ImportError:
        versions['cactus_commit'] = 'not_available'
    (out/'tool_provenance.json').write_text(json.dumps(versions, indent=2)+'\n', encoding='utf-8')
    write_table(out/'versions.tsv', [dict(process='PANGENOME_BUILD_QC', tool=tool, version=value) for tool,value in (
        ('cactus', versions['cactus_package']), ('cactus_commit', versions['cactus_commit']),
        ('toil', versions['toil_package']), ('verification_container', image))], ['process','tool','version'])
    result = dict(status='name_mapping_incomplete' if unknown else 'PASS', assemblies=len(records),
                  unknown_name_categories=len(unknown), full=full, clip=clip)
    (out/'clipping_audit.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    lines = ['# Graph representation and clipping', '',
        f"FULL support: {full['sequence_bp']:,} graph bp; CLIP biological graph: {clip['length']:,} graph bp.", '',
        'These are distinct graph-sequence totals, not assembly lengths or source-base retention percentages.', '',
        '| Assembly | Input FASTA bp | Retained path bp | Reported clip intervals, bp |',
        '|---|---:|---:|---:|']
    lines += [f"| {r['assembly']} | {r['source_fasta_bp']:,} | {r['retained_path_bp']:,} | {r['reported_clip_interval_bp']:,} |" for r in records]
    lines += ['', 'Retained path bp preserves repeated traversals and comes from Cactus sample-stats; '
        'it is different from the deduplicated sequence in sharing plots. Clip intervals are unioned within '
        'each tool-reported path name. They do not account for every sequence excluded before/during alignment, '
        'nor necessarily every later cleanup operation. The signed residual is an accounting diagnostic, '
        'not an estimate of deletion or clipping. Reference protection and assembly fragmentation can affect these values.', '',
        f"Tool-name mapping: {result['status']}. Unmapped tool names are explicit in the audit; reported clip totals are incomplete when names cannot be resolved.", '',
        '[Assembly accounting](assembly_clipping.tsv) · [Graph partitions](partition_retained_paths.tsv) · '
        '[FULL/CLIP totals](full_clip_graph.tsv) · [Unmapped names](unmapped_tool_names.tsv) · [Bundle inventory](clipping_bundle_inventory.tsv) · [Audit](clipping_audit.json)', '',
        '## Paths observed in FULL but not retained under the same sequence name in CLIP', '',
        'FULL W-record coordinate spans and CLIP retained path lengths are compared by biological assembly '
        'and source sequence name. FULL reference P-records are not assigned guessed lengths. '
        'This comparison can localize loss between graph stages; it does not certify a biological deletion.', '',
        '[Per-sequence FULL/CLIP evidence](full_clip_sequence_paths.tsv)', '',
        '| Assembly | Sequence | FULL W-coordinate span, bp | CLIP path bp |', '|---|---|---:|---:|']
    lines += [f"| {r['assembly']} | {r['sequence']} | {r['full_W_coordinate_span_bp']:,} | 0 |" for r in lost[:10]]
    if not lost:
        lines += ['No lost path names found among sequences with measured FULL W-coordinate spans.']
    lines += ['',
        '## Tool provenance', '',
        f"Verification image: `{image}`. Cactus package: `{versions['cactus_package']}`; commit: `{versions['cactus_commit']}`; Toil: `{versions['toil_package']}`.", '',
        'This separate check records versions in the verification image without rebuilding the graph. '
        'The original construction version record is preserved in [tool provenance](tool_provenance.json). '
        'The old `cactus --version` value uses Toil\'s version option and must not be labeled as the Cactus release.', '']
    (out/'clipping_report.md').write_text('\n'.join(lines), encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('bundle','full-gfa','clip-stats','identities','out','image','build-versions'):
        parser.add_argument('--'+key, required=True)
    args = parser.parse_args()
    summarize(args.bundle, args.full_gfa, args.clip_stats, args.identities, args.out, args.image, args.build_versions)


if __name__ == '__main__':
    main()
