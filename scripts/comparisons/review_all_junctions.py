#!/usr/bin/env python3
"""Inventory assessment boundaries; measure only explicitly selected suspects/controls."""
import argparse
from collections import defaultdict, Counter
import json
from pathlib import Path
import re
import shutil

from trace_chimera_origins import Fasta, discover, map_queries, run
from chimera_origin import sha, table, agp_rows
from assess_junction_batch import find_gaps, write_rows, verify_agp, hic_evidence
from junction_assessment import stable_id, validate_registry
from junction_focus import anchor_intervals
from junction_review import n_runs, compare_anchors, recommend, apply_reviews, select_evidence_rows


def specifications(a):
    if a.manifest:
        specs = json.loads(a.manifest.read_text())
    else:
        folder = a.results/'assembly/scaffold/yahs_round2'
        specs = []
        for path in sorted(folder.glob('*_round2_scaffolds.fa')):
            assembly = path.name.removesuffix('_round2_scaffolds.fa')
            match = re.fullmatch(r'(.+)_hap([12])', assembly)
            if not match:
                raise ValueError('Provide --manifest for nonstandard assembly identities')
            sample, hap = match.groups()
            specs.append(dict(assembly=assembly, sample=sample, haplotype=hap,
                              sister=sample+'_hap'+('2' if hap == '1' else '1'),
                              stages=discover(a.results, assembly, 'round2_only')))
    ids = [s['assembly'] for s in specs]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError('No assemblies or duplicate assembly IDs')
    by_id = {s['assembly']: s for s in specs}
    for s in specs:
        if sum(t['stage']=='assessment' for t in s['stages']) != 1:
            raise ValueError('Exactly one assessment stage is required')
        sister = s.get('sister', '')
        if sister:
            if sister == s['assembly'] or sister not in by_id or by_id[sister]['sample'] != s['sample']:
                raise ValueError('Sister must be a distinct explicitly paired assembly of the same individual')
    return specs


def inventory(a, specs):
    jobs, rows, audit = [], [], []
    seeds = validate_registry(table(a.registry)) if a.registry else []
    if set(r['assembly'] for r in seeds) - {s['assembly'] for s in specs}:
        raise ValueError('Registry includes an assembly outside the inventory')
    for i, spec in enumerate(specs):
        dest = a.out/f'assembly_{i}'; dest.mkdir()
        assessment = next(s for s in spec['stages'] if s['stage'] == 'assessment')
        digest = sha(Path(assessment['fasta']))
        fasta = Fasta(Path(assessment['fasta']), dest/'assessment.fa', a.samtools)
        print(f"Inventorying every scaffold in {spec['assembly']}", flush=True)
        gaps = find_gaps(spec['assembly'], fasta, spec['stages'], dest, a.samtools, set(fasta.lengths))
        write_rows(dest/'stage_inputs.tsv', [dict(stage=s['stage'], fasta=s['fasta'],
            agp=s.get('agp',''), agp_sha256=sha(Path(s['agp'])) if s.get('agp') and Path(s['agp']).is_file() else '.')
            for s in spec['stages']])
        by_interval = {(g['scaffold'], g['start'], g['end']): g for g in gaps}
        current = []
        for scaffold, length in fasta.lengths.items():
            sequence = fasta.fetch(scaffold, 0, length) if length else ''
            runs = n_runs(sequence)
            for lo, hi in runs:
                g = by_interval.pop((scaffold, lo, hi), None)
                current.append(dict(g or {}, id=stable_id(spec['assembly'], scaffold, lo, hi),
                    scaffold=scaffold, start=lo, end=hi,
                    kind='verified_agp_gap' if g else 'unresolved_N_run',
                    origin_stage=g['origin_stage'] if g else 'unknown',
                    left_component=g['left_component'] if g else '.',
                    right_component=g['right_component'] if g else '.'))
            audit.append(dict(assembly=spec['assembly'], scaffold=scaffold, length=length,
                              N_runs=len(runs), scanned=True))
        # AGP gaps can be adjacent to pre-existing Ns, producing a larger FASTA N run.
        for g in by_interval.values():
            current.append(dict(g, kind='verified_agp_gap'))
        if assessment.get('agp') and Path(assessment['agp']).is_file():
            agp = agp_rows(assessment['agp']); verify_agp(agp, fasta)
            for left, right in zip(agp, agp[1:]):
                if left['object'] == right['object'] and left['kind'] not in ('N','U') and right['kind'] not in ('N','U'):
                    pos = right['start']
                    current.append(dict(id=stable_id(spec['assembly'], right['object'], pos, pos),
                        scaffold=right['object'], start=pos, end=pos, kind='ungapped_agp_boundary',
                        origin_stage='assessment', left_component=left['component'], right_component=right['component']))
        for seed in (s for s in seeds if s['assembly'] == spec['assembly']):
            if seed['assessment_sha256'] != digest or seed['scaffold'] not in fasta.lengths or int(seed['end']) > fasta.lengths[seed['scaffold']]:
                raise ValueError('Registry differs from assessment sequence or coordinates')
            current.append(dict(seed, start=int(seed['start']), end=int(seed['end']),
                                kind='flagged_interval', origin_stage='requires_localization'))
        for r in current:
            r.update(assembly=spec['assembly'], sample=spec['sample'], haplotype=spec['haplotype'],
                     sister=spec.get('sister',''), assessment_sha256=digest, role='inventory',
                     hifi_chain_molecules='.', hic_status='not_measured',
                     decision='UNREVIEWED', reviewer='', reason='', evidence_for='', evidence_against='')
        rows += current
        jobs.append(dict(assembly=spec['assembly'], fasta=fasta, stages=spec['stages'],
                         dest=dest, selected=current, boundary_inventory=current, digest=digest, spec=spec))
    if len({r['id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate inventory/registry identifiers')
    write_rows(a.out/'scaffold_inventory.tsv', audit)
    return jobs, rows


def scope_evidence(a, jobs, rows):
    selections = [r for path in a.selection for r in table(path)]
    for packet in a.selection_packet:
        status = json.loads((packet/'status.json').read_text())
        if status.get('status') != 'SUCCESS':
            raise ValueError('Selection packet must be a completed assessment packet')
        for assembly in status['assemblies']:
            # Existing assessment tables have exact coordinates but carry the digest
            # in status.json rather than on each selected_gaps.tsv row.
            for source in table(packet/assembly['folder']/'selected_gaps.tsv'):
                if source['assembly'] != assembly['assembly']:
                    raise ValueError('Selection packet assembly mismatch')
                selections.append(dict(source, assessment_sha256=assembly['assessment_sha256']))
    if a.registry:
        selections.extend(dict(r, role='candidate_interval') for r in table(a.registry))
    chosen = select_evidence_rows(rows, selections)
    for job in jobs:
        job['selected'] = [r for r in chosen if r['assembly'] == job['assembly']]
        print(f"Evidence scope {job['assembly']}: {len(job['selected'])} selected / "
              f"{len(job['boundary_inventory'])} inventoried boundaries", flush=True)
    write_rows(a.out/'selected_junctions.tsv', chosen)
    return chosen


def reuse_packets(a, rows):
    lookup = {r['id']: r for r in rows}
    imported = []
    for packet in a.evidence_packet:
        status = json.loads((packet/'status.json').read_text())
        if status.get('status') != 'SUCCESS':
            raise ValueError('Evidence packet did not complete successfully')
        for assembly in status['assemblies']:
            folder = packet/assembly['folder']
            path = folder/'hifi_support.tsv'
            if not path.is_file():
                continue
            for source in table(path):
                r = lookup.get(source['id'])
                if not r:
                    continue
                if r['assembly'] != assembly['assembly'] or r['assessment_sha256'] != assembly['assessment_sha256'] or any(str(r[k]) != source[k] for k in ('scaffold','start','end')):
                    raise ValueError('Evidence packet coordinate/checksum mismatch')
                value = source.get('molecules_passing_chain_screen', '.') or '.'
                if r['hifi_chain_molecules'] != '.' and str(r['hifi_chain_molecules']) != value:
                    raise ValueError('Conflicting repeated measurements; select one evidence packet')
                r['hifi_chain_molecules'] = value
                imported.append(dict(id=r['id'], path=str(path.resolve()), sha256=sha(path)))
    write_rows(a.out/'reused_evidence.tsv', imported)


def sister_evidence(a, jobs):
    queries = a.out/'anchors.fa'; windows = []; mapping = {}; results = []
    with queries.open('w') as out:
        for job in jobs:
            gaps = job['selected']
            for r in gaps:
                for offset in a.offsets:
                    pair = []
                    for side, interval in anchor_intervals(r['start'], r['end'], job['fasta'].lengths[r['scaffold']], a.anchor_size, offset).items():
                        key = f"q{len(windows)}"
                        w = dict(query=key, id=r['id'], offset=offset, side=side, status='outside_scaffold')
                        if interval:
                            lo, hi = interval; seq = job['fasta'].fetch(r['scaffold'],lo,hi)
                            clo, chi = (lo,r['start']) if side=='left' else (r['end'],hi)
                            other = sum(g['id'] != r['id'] and g['kind'] != 'flagged_interval' and g['scaffold']==r['scaffold'] and g['start'] < chi and clo < g['end'] for g in job['boundary_inventory'])
                            w.update(start=lo,end=hi,status='clean' if set(seq)<=set('ACGT') and not other else 'ambiguous_or_other_boundary', other_boundaries=other)
                            out.write(f'>{key}\n{seq}\n')
                        windows.append(w); pair.append(w)
                    mapping[r['id'],offset] = pair
    write_rows(a.out/'anchor_windows.tsv', windows)
    for job in jobs:
        relevant = [r for j in jobs for r in j['selected'] if r['sister']==job['assembly']]
        if not relevant:
            continue
        print(f"Mapping sister anchors against {job['assembly']}", flush=True)
        wanted = {w['query'] for r in relevant for offset in a.offsets for w in mapping[r['id'],offset]}
        local = job['dest']/'sister_queries.fa'
        with queries.open() as src, local.open('w') as out:
            keep = False
            for line in src:
                if line.startswith('>'): keep = line[1:].strip() in wanted
                if keep: out.write(line)
        hits = map_queries(a, job['fasta'], local, job['dest']/'sister') if local.stat().st_size else []
        print(f"Sister mapping complete {job['assembly']}: {len(hits)} emitted alignments", flush=True)
        grouped = defaultdict(list)
        for h in hits: grouped[h['query']].append(h)
        for r in relevant:
            for offset in a.offsets:
                l, rr = mapping[r['id'],offset]
                results.append(dict(id=r['id'], sister=job['assembly'], offset=offset,
                    anchor_context='clean' if l['status']==rr['status']=='clean' else 'confounded',
                    expected_source_separation=r['end']-r['start']+2*offset,
                    **compare_anchors(grouped[l['query']],grouped[rr['query']])))
    write_rows(a.out/'sister_comparisons.tsv',results)
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--manifest',type=Path)
    p.add_argument('--registry',type=Path)
    p.add_argument('--selection',type=Path,action='append',default=[],
                   help='Explicit id/assembly/scaffold/start/end/assessment_sha256/role TSV; repeatable')
    p.add_argument('--selection-packet',type=Path,action='append',default=[],
                   help='Completed junction-assessment packet supplying bounded selected gaps and controls')
    p.add_argument('--evidence-packet',type=Path,action='append',default=[])
    p.add_argument('--reviews',type=Path)
    p.add_argument('--inventory-only',action='store_true')
    p.add_argument('--hic',action='store_true')
    p.add_argument('--hic-inputs',type=Path)
    p.add_argument('--contact-flank',type=int,default=100000)
    p.add_argument('--anchor-size',type=int,default=10000)
    p.add_argument('--offsets',type=int,nargs='+',default=[0,25000,100000])
    p.add_argument('--threads',type=int,default=8)
    p.add_argument('--index-bases',type=int,default=16000000000)
    p.add_argument('--samtools',default='samtools')
    p.add_argument('--minimap2',default='minimap2')
    a=p.parse_args()
    if not a.inventory_only and not (a.selection or a.selection_packet or a.registry):
        p.error('Evidence analysis requires --selection, --selection-packet or --registry; '
                'use --inventory-only for unrestricted inventory')
    if min(a.anchor_size,a.contact_flank,a.threads,a.index_bases)<=0 or min(a.offsets)<0 or len(set(a.offsets))!=len(a.offsets):
        p.error('Positive sizes/threads and distinct nonnegative offsets required')
    a.out.mkdir(parents=True,exist_ok=False)
    specs=specifications(a)
    (a.out/'input_manifest.json').write_text(json.dumps(specs,indent=2))
    (a.out/'parameters.json').write_text(json.dumps(vars(a),default=str,indent=2))
    run([a.samtools,'--version'],stdout=a.out/'samtools.version.txt')
    if not a.inventory_only:
        run([a.minimap2,'--version'],stdout=a.out/'minimap2.version.txt')
    a.hic_input_map={r['assembly']:r for r in table(a.hic_inputs)} if a.hic_inputs else {}
    jobs,rows=inventory(a,specs)
    write_rows(a.out/'junction_inventory.tsv',rows)
    if not a.inventory_only:
        rows=scope_evidence(a,jobs,rows)
    reuse_packets(a,rows)
    sisters=[] if a.inventory_only else sister_evidence(a,jobs)
    by_id=defaultdict(list)
    for s in sisters: by_id[s['id']].append(s)
    hic=[]
    for job in jobs:
        status=hic_evidence(a,job,job['dest']) if a.hic and not a.inventory_only and job['selected'] else dict(status='not_requested')
        hic.append(dict(assembly=job['assembly'],**status))
        for r in job['selected']: r['hic_status']=status['status']
    for r in rows:
        r['sister_states']=';'.join(sorted({s['status'] for s in by_id[r['id']]})) or 'not_assessed'
        r['recommendation'],r['recommendation_reason'],r['next_informative_check']=recommend(r,by_id[r['id']])
    if a.reviews: apply_reviews(rows,table(a.reviews))
    write_rows(a.out/'join_decisions.tsv',rows,list(dict.fromkeys(k for r in rows for k in r)))
    for path in (Path(__file__),Path(__file__).parents[2]/'py_scripts/junction_review.py'):
        shutil.copy2(path,a.out/path.name)
    (a.out/'status.json').write_text(json.dumps(dict(status='SUCCESS',cut_authorized=False,
        junctions=len(rows),assemblies=[dict(assembly=j['assembly'],sha256=j['digest']) for j in jobs],
        recommendations=dict(Counter(r['recommendation'] for r in rows)),hic=hic),indent=2))


if __name__=='__main__': main()
