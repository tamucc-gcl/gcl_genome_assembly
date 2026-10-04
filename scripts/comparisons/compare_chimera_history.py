#!/usr/bin/env python3
"""Compare assessment assemblies to historical assemblies without changing either."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import shutil

from trace_chimera_origins import Fasta, discover, map_queries, run
from chimera_origin import exact_projection, paf_rows, sha, table, write_table


def union_length(intervals):
    end, total = -1, 0
    for lo, hi in sorted(intervals):
        if lo < 0 or hi <= lo:
            raise ValueError('Invalid interval')
        total += max(0, hi-max(lo, end))
        end = max(end, hi)
    return total


def summarize(hits):
    groups = defaultdict(list)
    for h in hits:
        groups[h['query'], h['target'], h['strand']].append(h)
    return [dict(query=q, target=t, strand=s, blocks=len(hs),
                 query_covered_bp=union_length((h['query_start'], h['query_end']) for h in hs),
                 target_covered_bp=union_length((h['target_start'], h['target_end']) for h in hs),
                 query_length=hs[0]['query_length'], target_length=hs[0]['target_length'])
            for (q, t, s), hs in sorted(groups.items())]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--old', action='append', required=True, help='ID=FASTA; repeat for all historical haplotypes')
    p.add_argument('--assembly', action='append', required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--flank', type=int, default=100000)
    p.add_argument('--anchor', type=int, default=1000)
    p.add_argument('--index-bases', type=int, default=8000000000)
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--minimap2', default='minimap2')
    a = p.parse_args()
    if min(a.threads, a.flank, a.anchor, a.index_bases) <= 0 or a.anchor > a.flank:
        p.error('Positive sizes required, with anchor <= flank')
    if len(set(a.assembly)) != len(a.assembly):
        p.error('Duplicate assembly IDs')
    old_specs = [s.split('=', 1) for s in a.old]
    if any(len(s) != 2 or not s[0] for s in old_specs) or len({s[0] for s in old_specs}) != len(old_specs):
        p.error('Historical assemblies require unique ID=FASTA entries')
    a.out.mkdir(parents=True, exist_ok=False)
    for src in (Path(__file__), Path(__file__).with_name('trace_chimera_origins.py'),
                Path(__file__).resolve().parents[2]/'py_scripts/chimera_origin.py'):
        shutil.copy2(src, a.out/src.name)
    provenance = dict(arguments=vars(a).copy(), inputs=[], coordinate_system='zero_based_half_open')
    provenance['arguments'] = {k: str(v) if isinstance(v, Path) else v for k, v in provenance['arguments'].items()}
    for tool, flag, label in ((a.minimap2, '--version', 'minimap2'), (a.samtools, '--version', 'samtools')):
        run([tool, flag], stdout=a.out/f'{label}.version.txt')
    old = []
    for n, (label, filename) in enumerate(old_specs):
        path = Path(filename).resolve(strict=True)
        fasta = Fasta(path, a.out/f'old_{n}.fa', a.samtools)
        old.append((label, fasta))
        provenance['inputs'].append(dict(role='historical', id=label, path=str(path), sha256=sha(path)))
    for n, assembly in enumerate(a.assembly):
        if Path(assembly).name != assembly or '/' in assembly or '\\' in assembly:
            raise ValueError('Assembly ID must be a filename component')
        dest = a.out/f'current_{n}'
        dest.mkdir()
        calls_path = a.results/'assembly/chimeras'/f'{assembly}.chimeric_joins.tsv'
        calls = [r for r in table(calls_path) if r.get('chromosome_member') == 'yes' and r.get('transition_lo', '.') != '.']
        if not calls or any(r['assembly'] != assembly for r in calls):
            raise ValueError('Expected chromosome candidates for requested assembly')
        scopes = {r['join_scope'] for r in calls}
        if len(scopes) != 1 or any(r['coordinate_stage'] != 'pre_finishing' for r in calls):
            raise ValueError('Mixed scaffold scopes')
        specs = discover(a.results, assembly, scopes.pop())
        path = Path(next(s['fasta'] for s in specs if s['stage'] == 'assessment')).resolve(strict=True)
        digest = sha(path)
        if {r['assembly_sha256'] for r in calls} != {digest}:
            raise ValueError('Current assessment checksum mismatch')
        current = Fasta(path, dest/'assessment.fa', a.samtools)
        provenance['inputs'].append(dict(role='current_assessment', id=assembly, path=str(path), sha256=digest))
        shutil.copy2(calls_path, dest/'candidates.tsv')
        windows, intervals = dest/'windows.fa', []
        with windows.open('w') as out:
            for i, r in enumerate(calls, 1):
                lo, hi = int(r['transition_lo']), int(r['transition_hi'])
                scaffold = r['scaffold']
                start, end = max(0, lo-a.flank), min(current.lengths[scaffold], hi+a.flank)
                if not start <= lo < hi <= end:
                    raise ValueError('Invalid candidate interval')
                label = f'interval_{i:03d}'
                out.write(f'>{label}\n{current.fetch(scaffold, start, end)}\n')
                intervals.append(dict(candidate=label, scaffold=scaffold, start=start, end=end, lo=lo, hi=hi))
        write_table(dest/'intervals.tsv', ['candidate', 'scaffold', 'start', 'end', 'lo', 'hi'], intervals)
        for j, (old_id, target) in enumerate(old):
            print(f'{assembly} against {old_id}', flush=True)
            prefix = dest/f'old_{j}'
            # Full assembly context includes secondary matches; no base-level CIGAR is needed here.
            if sum(target.lengths.values()) >= a.index_bases:
                raise ValueError('Historical assembly exceeds single-index resource budget')
            cmd = [a.minimap2, '-x', 'asm5', '-I', str(a.index_bases), '--secondary=yes',
                   '-N', '50', '-p', '0.5', '-t', str(a.threads), str(target.path), str(current.path)]
            Path(str(prefix)+'.whole.command.json').write_text(json.dumps(cmd))
            run(cmd, stdout=Path(str(prefix)+'.whole.paf'), stderr=Path(str(prefix)+'.whole.log'))
            hits = paf_rows(Path(str(prefix)+'.whole.paf'))
            write_table(Path(str(prefix)+'.chromosome_matches.tsv'),
                        ['query', 'target', 'strand', 'blocks', 'query_covered_bp', 'target_covered_bp', 'query_length', 'target_length'], summarize(hits))
            window_hits = map_queries(a, target, windows, dest/f'old_{j}_windows')
            anchors = []
            for r in intervals:
                for side, lo, hi in [('left', r['lo']-a.anchor, r['lo']), ('right', r['hi'], r['hi']+a.anchor)]:
                    found = False
                    if r['start'] <= lo < hi <= r['end'] and set(current.fetch(r['scaffold'], lo, hi)) <= set('ACGT'):
                        for h in window_hits:
                            if h['query'] != r['candidate']:
                                continue
                            projection = exact_projection(h, lo-r['start'], hi-r['start'])
                            if projection:
                                found = True
                                anchors.append(dict(candidate=r['candidate'], side=side, old_id=old_id,
                                                    target=h['target'], start=projection[0], end=projection[1],
                                                    strand=h['strand'], status='exact_anchor_in_reported_alignment'))
                    if not found:
                        anchors.append(dict(candidate=r['candidate'], side=side, old_id=old_id, target='.', start='.', end='.', strand='.', status='not_established'))
            write_table(Path(str(prefix)+'.anchors.tsv'), ['candidate', 'side', 'old_id', 'target', 'start', 'end', 'strand', 'status'], anchors)
    (a.out/'provenance.json').write_text(json.dumps(provenance, indent=2))
    (a.out/'README.txt').write_text('Historical assemblies are comparison evidence, not truth. All current/old combinations are aligned; haplotype labels need not correspond. Whole PAFs preserve secondary matches. Coverage is a union per chromosome pair/orientation, not unique coverage; do not sum across partners. Windows PAFs contain base-resolved alignments. Missing exact anchors may reflect divergence, edits or alignment limitations; they do not establish absence. Same-chromosome anchors do not alone prove intervening continuity. No cut decisions are made.\n')
    (a.out/'status.json').write_text(json.dumps(dict(status='SUCCESS', cut_authorized=False)))


if __name__ == '__main__':
    main()
