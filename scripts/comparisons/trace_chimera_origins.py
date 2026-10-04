#!/usr/bin/env python3
"""Trace candidate intervals through retained stages; never edit assemblies or emit cuts."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'py_scripts'))
from chimera_origin import (agp_rows, component_interval, exact_projection, exact_span, paf_rows, read_support,
                            sha, source_to_query, table, write_table)


def run(command, stdout=None, stderr=None):
    with open(stderr or os.devnull, 'w') as err:
        if stdout:
            with open(stdout, 'w') as out:
                subprocess.run([str(v) for v in command], stdout=out, stderr=err, check=True)
        else:
            subprocess.run([str(v) for v in command], stderr=err, check=True)


class Fasta:
    def __init__(self, path, local, samtools):
        self.path, self.samtools = local, samtools
        local.symlink_to(path.resolve())
        index_log = local.with_suffix('.index.log')
        run([samtools, 'faidx', local], stderr=index_log)
        if 'duplicate sequence' in index_log.read_text().lower():
            raise ValueError('Duplicate FASTA identifier reported by samtools')
        self.lengths = {}
        for line in Path(str(local)+'.fai').read_text().splitlines():
            f = line.split('\t')
            if f[0] in self.lengths:
                raise ValueError('Duplicate FASTA identifier')
            self.lengths[f[0]] = int(f[1])

    def fetch(self, name, lo, hi):
        if not 0 <= lo < hi <= self.lengths[name]:
            raise ValueError(f'Invalid FASTA interval: {name}:{lo}-{hi}')
        value = subprocess.check_output([self.samtools, 'faidx', str(self.path),
                                         f'{name}:{lo+1}-{hi}'], text=True)
        sequence = ''.join(value.splitlines()[1:]).upper()
        if len(sequence) != hi-lo:
            raise ValueError('FASTA extraction length mismatch')
        return sequence


def discover(results, assembly, scope):
    """Layout adapter only. The tracing algorithm takes an explicit stage manifest."""
    match = re.fullmatch(r'(.+)_(hap[12]|primary)', assembly)
    if not match:
        raise ValueError('Use --manifest for assembly IDs outside the published layout convention')
    sample, hap = match.groups()
    contig = results/'assembly/contig'
    scaffold = results/'assembly/scaffold'
    choices = [
        ('raw_contigs', contig/'hifiasm'/f'{sample}.{hap}.p_ctg.fasta', None),
        ('organelle_filtered', contig/'organelle_filtered'/f'{assembly}.organelle_filtered.fasta', None),
        ('purged', contig/'purge_dups'/f'{assembly}.purged.fa', None),
        ('contig_corrected', contig/'misassembly_correction'/f'{assembly}_corrected.fasta', None),
        ('round1_input', contig/'decontam'/f'{assembly}.decontaminated.fasta', None),
        ('round1', scaffold/'yahs'/f'{assembly}_scaffolds.fa',
         scaffold/'yahs'/f'{assembly}_scaffolds_final.agp'),
    ]
    if scope == 'round2_only':
        choices.extend([
            ('scaffold_corrected', scaffold/'misassembly_correction'/f'{assembly}_corrected.fasta', None),
            ('assessment', scaffold/'yahs_round2'/f'{assembly}_round2_scaffolds.fa',
             scaffold/'yahs_round2'/f'{assembly}_round2_scaffolds_final.agp')])
    elif scope == 'round1_only':
        choices[-1] = ('assessment', choices[-1][1], choices[-1][2])
    else:
        raise ValueError('Unknown assessment round: use an explicit manifest')
    return [dict(stage=name, fasta=str(fa.resolve()), agp=str(agp.resolve()) if agp else '',
                 graph=str((contig/'hifiasm'/f'{sample}.{hap}.p_ctg.gfa').resolve())
                 if name == 'raw_contigs' else '') for name, fa, agp in choices]


def map_queries(args, target, query, prefix):
    if sum(target.lengths.values()) >= args.index_bases:
        raise ValueError('Target exceeds configured single-index budget; increase --index-bases/resources')
    command = [args.minimap2, '-x', 'asm5', '-I', str(args.index_bases), '-c', '--eqx',
               '--secondary=yes', '-N', '50', '-p', '0.5', '-t', str(args.threads), target.path, query]
    run(command, prefix.with_suffix('.paf'), prefix.with_suffix('.log'))
    prefix.with_suffix('.command.json').write_text(json.dumps([str(v) for v in command], indent=2)+'\n')
    return paf_rows(prefix.with_suffix('.paf'))


def cached_inventory(trace, work_root, out):
    """Inspect just task roots named in the trace, not the whole work filesystem."""
    rows = []
    for task in table(trace):
        if not any(x in task['name'] for x in (':HIFIASM (', ':SCAFFOLD_HIC')):
            continue
        match = re.fullmatch(r'([0-9a-f]{2})/([0-9a-f]+)', task['hash'])
        if not match:
            raise ValueError('Invalid trace task hash')
        parent = work_root/match[1]
        roots = list(parent.glob(match[2]+'*'))
        roots = [p for p in roots if p.is_dir()]
        if len(roots) != 1:
            rows.append(dict(task=task['name'], path=str(parent/match[2]),
                             status='missing_or_ambiguous_task_directory', bytes='.'))
            continue
        found = []
        for path in roots[0].iterdir():
            if path.suffix in ('.gfa', '.bin', '.bed', '.agp'):
                found.append(path)
                rows.append(dict(task=task['name'], path=str(path.resolve()),
                                 status='retained' if path.is_file() else 'broken_link',
                                 bytes=path.stat().st_size if path.is_file() else '.'))
        if not found:
            rows.append(dict(task=task['name'], path=str(roots[0]), status='no_retained_artifacts', bytes='.'))
    write_table(out, ['task', 'path', 'status', 'bytes'], rows)


def graph_summary(path, targets, out):
    """Read topology metadata only; never put gigabase segment sequences in reports."""
    counts, seen, links = Counter(), set(), []
    with open(path) as handle:
        for line in handle:
            # S records can contain whole chromosomes: discard sequence immediately.
            f = line.rstrip().split('\t')
            counts[f[0]] += 1
            if f[0] == 'S' and len(f) >= 3 and f[1] in targets:
                seen.add(f[1])
            if f[0] == 'L' and len(f) >= 6 and (f[1] in targets or f[3] in targets):
                links.append(dict(source=f[1], source_orientation=f[2], target=f[3],
                                  target_orientation=f[4], overlap=f[5]))
    write_table(out/'raw_graph_links.tsv', ['source', 'source_orientation', 'target',
                                           'target_orientation', 'overlap'], links)
    (out/'raw_graph_summary.json').write_text(json.dumps(dict(path=str(path), records=counts,
        matched_segments=sorted(seen), requested_segments=sorted(targets),
        interpretation='Contig graph topology; missing links do not establish absent biological adjacency. '
                       'Unitig/read paths are not reconstructed here.'), indent=2)+'\n')


def partner_summary(path, intervals, lengths, flank, out):
    """Raw contacts by partner scaffold for each interval's external flanks.

    No distance/bias normalization or library independence is implied. Aggregate
    pairs lack library labels; explicit library inputs are needed for that test.
    """
    counts, by_scaffold = Counter(), defaultdict(list)
    for key, iv in intervals.items():
        by_scaffold[iv['scaffold']].append((key, 'left', max(0, iv['lo']-flank), iv['lo']))
        by_scaffold[iv['scaffold']].append((key, 'right', iv['hi'],
                                           min(lengths[iv['scaffold']], iv['hi']+flank)))
    read = 0
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip().split('\t')
            if len(f) < 8:
                raise ValueError('Malformed alternative partner pairs')
            read += 1
            if f[7] != 'UU' or f[1] == f[3]:
                continue
            for chrom, pos, partner in ((f[1], int(f[2])-1, f[3]), (f[3], int(f[4])-1, f[1])):
                for key, side, lo, hi in by_scaffold.get(chrom, []):
                    if lo <= pos < hi:
                        counts[key, side, partner] += 1
    rows = []
    for chrom, windows in by_scaffold.items():
        for key, side, lo, hi in windows:
            partners = [(partner, count) for (k, s, partner), count in counts.items() if k == key and s == side]
            for partner, count in sorted(partners, key=lambda x: (-x[1], x[0])) or [('.', 0)]:
                rows.append(dict(candidate=key, side=side, scaffold=chrom, start=lo, end=hi,
                                 partner=partner, contacts=count, library='pooled_unknown'))
    write_table(out, ['candidate', 'side', 'scaffold', 'start', 'end', 'partner', 'contacts', 'library'], rows)
    return read


def process_assembly(args, assembly, calls, stages, dest):
    dest.mkdir()
    (dest/'requested_inputs.json').write_text(json.dumps(dict(assembly=assembly, stages=stages,
                                                           calls=calls), indent=2)+'\n')
    digests = {r['assembly_sha256'] for r in calls}
    if len(digests) != 1 or not re.fullmatch('[0-9a-f]{64}', next(iter(digests))):
        raise ValueError('Calls must identify exactly one assessment SHA256')
    digest = next(iter(digests))
    if len({s['stage'] for s in stages}) != len(stages) or sum(s['stage'] == 'assessment' for s in stages) != 1:
        raise ValueError('Manifest needs unique stage names and exactly one assessment')
    assessment_spec = next(s for s in stages if s['stage'] == 'assessment')
    if sha(assessment_spec['fasta']) != digest:
        raise ValueError(f'{assembly}: assessment checksum mismatch')
    assessment = Fasta(Path(assessment_spec['fasta']), dest/'assessment.fa', args.samtools)
    intervals, sequences = {}, {}
    for i, row in enumerate(calls, 1):
        key = f'interval_{i:03d}'
        lo, hi, scaffold = int(row['transition_lo']), int(row['transition_hi']), row['scaffold']
        if not 0 <= lo <= hi <= assessment.lengths[scaffold]:
            raise ValueError('Candidate bounds outside assessment')
        start, end = max(0, lo-args.flank), min(assessment.lengths[scaffold], hi+args.flank)
        intervals[key] = dict(scaffold=scaffold, lo=lo, hi=hi, start=start, end=end,
                              left_chrom=row.get('left_chrom', '.'), right_chrom=row.get('right_chrom', '.'))
        sequences[key] = assessment.fetch(scaffold, start, end)
    queries = dest/'intervals.fa'
    queries.write_text(''.join(f'>{key}\n{seq}\n' for key, seq in sequences.items()))
    ledger, boundaries, components, stage_audit, flank_placements, raw_targets = [], [], [], [], [], set()
    for n, spec in enumerate(stages):
        stage = spec['stage']
        if not re.fullmatch('[A-Za-z0-9_-]+', stage):
            raise ValueError('Unsafe stage name')
        path = Path(spec['fasta'])
        if not path.is_file():
            stage_audit.append(dict(stage=stage, path=str(path), sha256='.', status='missing',
                                    agp_path=spec.get('agp', ''), agp_sha256='.', agp_status='not_assessed'))
            for key in intervals:
                ledger.append(dict(candidate=key, stage=stage, status='stage_missing', exact_placements=0,
                                   target='.', target_start='.', target_end='.', strand='.', mapq='.',
                                   includes_N='.', action='REVIEW'))
            continue
        print(f'{assembly}: tracing {stage}', flush=True)
        target = assessment if stage == 'assessment' else Fasta(path, dest/f'stage_{n}.fa', args.samtools)
        agp_available = bool(spec.get('agp') and Path(spec['agp']).is_file())
        stage_audit.append(dict(stage=stage, path=str(path), sha256=digest if stage == 'assessment' else sha(path),
            status='available', agp_path=spec.get('agp', ''),
            agp_sha256=sha(spec['agp']) if agp_available else '.',
            agp_status='available' if agp_available else 'missing' if spec.get('agp') else 'not_applicable'))
        agp = agp_rows(spec['agp']) if agp_available else []
        if agp_available:
            lengths = {}
            for row in agp:
                lengths[row['object']] = row['end']
            if lengths != target.lengths:
                raise ValueError(f'{assembly} {stage}: AGP and FASTA object lengths disagree')
        if stage == 'assessment':
            hits = [dict(query=k, query_length=len(sequences[k]), query_start=0, query_end=len(sequences[k]),
                         target=iv['scaffold'], target_length=assessment.lengths[iv['scaffold']],
                         target_start=iv['start'], target_end=iv['end'], strand='+', mapq=60,
                         cigar=f'{len(sequences[k])}=', alignment_type='identity') for k, iv in intervals.items()]
        else:
            hits = map_queries(args, target, queries, dest/f'{stage}_alignment')
        by_query = defaultdict(list)
        for hit in hits:
            if hit['query'] not in sequences:
                raise ValueError('Unknown query in PAF')
            by_query[hit['query']].append(hit)
        for key, iv in intervals.items():
            for side, alo, ahi in (('left', iv['lo']-args.anchor, iv['lo']),
                                   ('right', iv['hi'], iv['hi']+args.anchor)):
                placements = {}
                if (iv['start'] <= alo < ahi <= iv['end'] and not
                        re.search('[^ACGT]', sequences[key][alo-iv['start']:ahi-iv['start']])):
                    for h in by_query[key]:
                        projection = exact_projection(h, alo-iv['start'], ahi-iv['start'])
                        if projection:
                            placements[(h['target'],) + projection + (h['strand'],)] = h['mapq']
                for placement, quality in placements.items():
                    flank_placements.append(dict(candidate=key, stage=stage, side=side,
                        assessment_start=alo, assessment_end=ahi, target=placement[0],
                        target_start=placement[1], target_end=placement[2], strand=placement[3],
                        mapq=quality, reported_exact_placements=len(placements)))
                if not placements:
                    flank_placements.append(dict(candidate=key, stage=stage, side=side,
                        assessment_start=alo, assessment_end=ahi, target='.', target_start='.',
                        target_end='.', strand='.', mapq='.', reported_exact_placements=0))
            # Deduplicate identical emitted placements, not overlapping alternative hits.
            exact = {}
            for h in by_query[key]:
                if exact_span(h, sequences[key], target.fetch):
                    exact[h['target'], h['target_start'], h['target_end'], h['strand']] = h
            status = ('exact_sequence_present' if len(exact) == 1 else
                      'multiple_exact_placements' if exact else 'no_exact_full_window_placement')
            for h in list(exact.values()) or [None]:
                ledger.append(dict(candidate=key, stage=stage, status=status, exact_placements=len(exact),
                    target=h['target'] if h else '.', target_start=h['target_start'] if h else '.',
                    target_end=h['target_end'] if h else '.', strand=h['strand'] if h else '.',
                    mapq=h['mapq'] if h else '.', includes_N='yes' if 'N' in sequences[key] else 'no', action='REVIEW'))
                if h is None:
                    continue
                if stage == 'raw_contigs':
                    raw_targets.add(h['target'])
                for component in agp:
                    if component['object'] != h['target']:
                        continue
                    overlap_lo = max(component['start'], h['target_start'])
                    overlap_hi = min(component['end'], h['target_end'])
                    if component['kind'] not in ('N', 'U') and overlap_lo < overlap_hi:
                        clo, chi = component_interval(component, overlap_lo, overlap_hi)
                        qlo, qhi = source_to_query(h, overlap_lo, overlap_hi)
                        components.append(dict(candidate=key, stage=stage, object=component['object'],
                            component=component['component'], component_start=clo, component_end=chi,
                            assessment_scaffold=iv['scaffold'], assessment_start=qlo+iv['start'],
                            assessment_end=qhi+iv['start'],
                            assessment_to_component_orientation='+' if h['strand'] == component['orientation'] else '-',
                            placement_ambiguity='yes' if len(exact) > 1 else 'not_observed'))
                    # Emit whole gaps or component-end points only within the verified window.
                    spans = [(component['start'], component['end'], 'agp_gap')] if component['kind'] in ('N', 'U') else [
                        (component['start'], component['start'], 'component_start'),
                        (component['end'], component['end'], 'component_end')]
                    for lo, hi, kind in spans:
                        if not h['target_start'] <= lo <= hi <= h['target_end']:
                            continue
                        if kind == 'agp_gap' and set(target.fetch(h['target'], lo, hi)) != {'N'}:
                            raise ValueError('AGP gap disagrees with the verified stage sequence')
                        qlo, qhi = source_to_query(h, lo, hi)
                        alo, ahi = qlo+iv['start'], qhi+iv['start']
                        boundaries.append(dict(candidate=key, stage=stage, kind=kind,
                            source_object=component['object'], source_start=lo, source_end=hi,
                            component=component['component'], source_orientation=component['orientation'],
                            assessment_scaffold=iv['scaffold'], assessment_start=alo, assessment_end=ahi,
                            inside_transition='yes' if alo <= iv['hi'] and ahi >= iv['lo'] else 'no',
                            coordinate_status='literal_full_window_verified',
                            placement_ambiguity='yes' if len(exact) > 1 else 'not_observed', action='REVIEW'))
            # Detailed alignments are retained even when no exact full window exists.
            # Do not infer first creation from absent/changed/multimapping alignments.
    write_table(dest/'origin_ledger.tsv', ['candidate', 'stage', 'status', 'exact_placements', 'target',
        'target_start', 'target_end', 'strand', 'mapq', 'includes_N', 'action'], ledger)
    write_table(dest/'stage_inputs.tsv', ['stage', 'path', 'sha256', 'status', 'agp_path', 'agp_sha256', 'agp_status'], stage_audit)
    write_table(dest/'flank_placements.tsv', ['candidate', 'stage', 'side', 'assessment_start', 'assessment_end',
        'target', 'target_start', 'target_end', 'strand', 'mapq', 'reported_exact_placements'], flank_placements)
    write_table(dest/'source_boundaries.tsv', ['candidate', 'stage', 'kind', 'source_object',
        'source_start', 'source_end', 'component', 'source_orientation', 'assessment_scaffold',
        'assessment_start', 'assessment_end', 'inside_transition', 'coordinate_status', 'placement_ambiguity', 'action'], boundaries)
    write_table(dest/'source_components.tsv', ['candidate', 'stage', 'object', 'component',
        'component_start', 'component_end', 'assessment_scaffold', 'assessment_start', 'assessment_end',
        'assessment_to_component_orientation', 'placement_ambiguity'], components)
    graph = next((s.get('graph') for s in stages if s['stage'] == 'raw_contigs'), None)
    if graph and Path(graph).is_file():
        graph_summary(Path(graph), raw_targets, dest)
    evidence_checks(args, assembly, digest, assessment, intervals, sequences, boundaries, dest)
    pair_file = args.results/'assembly/chimeras/evidence'/f'{assembly}.alternative_partners.pairs.gz'
    contacts = partner_summary(pair_file, intervals, assessment.lengths, args.contact_flank,
                               dest/'alternative_partner_counts.tsv') if pair_file.is_file() else None
    (dest/'provenance.json').write_text(json.dumps(dict(assembly=assembly, assessment_sha256=digest,
        coordinate_stage='pre_finishing', intervals=intervals, stages=stages,
        parameters=vars_for_json(args), alternative_pairs_read=contacts,
        decision='REVIEW', cut_authorized=False), indent=2)+'\n')
    return dict(assembly=assembly, intervals=len(intervals), assessment_sha256=digest,
                earliest_exact_stage={key: next((r['stage'] for r in ledger if r['candidate'] == key and
                    r['exact_placements'] > 0), '.') for key in intervals})


def vars_for_json(args):
    return {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}


def evidence_checks(args, assembly, digest, assessment, intervals, sequences, boundaries, dest):
    """Screen whole uncertainty intervals and real boundaries; never invent midpoint cuts."""
    hypotheses = []
    for key, iv in intervals.items():
        positions = {(iv['lo'], iv['hi'], 'uncertainty_interval')}
        # Assessment N-runs remain a distinct, explicitly sequence-derived hypothesis.
        for match in re.finditer('N+', sequences[key]):
            lo, hi = match.start()+iv['start'], match.end()+iv['start']
            if lo <= iv['hi'] and hi >= iv['lo']:
                positions.add((lo, hi, 'assessment_N_gap'))
        for row in boundaries:
            if row['candidate'] == key and row['inside_transition'] == 'yes' and row['placement_ambiguity'] == 'not_observed':
                positions.add((row['assessment_start'], row['assessment_end'], 'verified_source_boundary'))
        for lo, hi, kind in sorted(positions):
            left, right = lo-args.anchor, hi+args.anchor
            status = 'available' if 0 <= left < lo <= hi < right <= assessment.lengths[iv['scaffold']] else 'short_flank'
            lseq = assessment.fetch(iv['scaffold'], left, lo) if status == 'available' else ''
            rseq = assessment.fetch(iv['scaffold'], hi, right) if status == 'available' else ''
            if status == 'available' and re.search('[^ACGT]', lseq+rseq):
                status = 'ambiguous_anchor_bases'
            hypotheses.append(dict(id=f'H{len(hypotheses)+1:05d}', candidate=key, kind=kind,
                scaffold=iv['scaffold'], start=lo, end=hi, left_anchor_start=left,
                right_anchor_end=right, anchor_status=status, left_sequence=lseq,
                right_sequence=rseq, cut_authorized='no'))
    anchors = dest/'anchors.fa'
    anchors.write_text(''.join(f'>{h["id"]}_{side}\n{h[side+"_sequence"]}\n' for h in hypotheses
                              if h['anchor_status'] == 'available' for side in ('left', 'right')))
    hits = map_queries(args, assessment, anchors, dest/'anchor_specificity') if anchors.stat().st_size else []
    placements, intended, emitted = defaultdict(set), set(), Counter()
    by_anchor = {h['id']+'_'+side: (h['scaffold'],
                 h['left_anchor_start'] if side == 'left' else h['end'],
                 h['start'] if side == 'left' else h['right_anchor_end'])
                 for h in hypotheses for side in ('left', 'right')}
    for hit in hits:
        emitted[hit['query']] += 1
        if hit['query_end']-hit['query_start'] >= 0.95*hit['query_length'] and hit['matches']/hit['block'] >= 0.98:
            placements[hit['query']].add((hit['target'], hit['target_start'], hit['target_end'], hit['strand']))
        if (hit['strand'] == '+' and hit['query_start'] == 0 and hit['query_end'] == args.anchor and
                hit['matches'] == hit['block'] == args.anchor and
                (hit['target'], hit['target_start'], hit['target_end']) == by_anchor[hit['query']]):
            intended.add(hit['query'])
    context = None
    if args.hifi_results:
        path = args.hifi_results/'assembly/chimeras/sequence_context'/f'{assembly}.sequence_context/provenance.json'
        if path.is_file():
            meta = json.loads(path.read_text())
            if meta.get('hifi'):
                if (meta.get('assembly') != assembly or meta.get('sha256') != digest or
                        meta.get('coordinate_stage') != 'pre_finishing'):
                    raise ValueError('Retained HiFi context has a different assessment coordinate frame')
                context = (path.parent, meta)
                shutil.copy2(path, dest/'retained_hifi_provenance.json')
    for h in hypotheses:
        for side in ('left', 'right'):
            reported = placements[h['id']+'_'+side]
            h[side+'_reported_placements'] = len(reported)
        h['intended_anchors_reported'] = 'yes' if all(h['id']+'_'+side in intended for side in ('left', 'right')) else 'no'
        h['possible_secondary_saturation'] = 'yes' if any(emitted[h['id']+'_'+side] >= 51 for side in ('left', 'right')) else 'no'
        h['anchor_specificity'] = ('multiple_reported_placements' if max(h['left_reported_placements'], h['right_reported_placements']) > 1 else
                                  'single_reported_placement_each' if min(h['left_reported_placements'], h['right_reported_placements']) == 1
                                  and h['intended_anchors_reported'] == 'yes' and h['possible_secondary_saturation'] == 'no' else
                                  'not_established')
        h['hifi_status'] = 'not_available'
        if h['anchor_status'] != 'available':
            h['hifi_status'] = 'anchors_unassessable'
        elif context:
            root, meta = context
            h['hifi_status'] = 'no_retained_interval_covering_anchors'
            matches = [name for name, iv in meta['intervals'].items() if iv['scaffold'] == h['scaffold'] and
                       iv['start'] <= h['left_anchor_start'] and iv['end'] >= h['right_anchor_end']]
            if matches:
                path = root/(matches[0]+'.sam')
                h['hifi_status'] = 'retained_SAM_missing'
                if path.is_file():
                    with path.open() as handle:
                        h.update(read_support(handle, h['scaffold'], h['left_anchor_start'], h['right_anchor_end']))
                    h['hifi_status'] = 'retained_regional_SAM_screened'
                    h['hifi_source'] = str(path.resolve())
                    h['hifi_source_sha256'] = sha(path)
    fields = ['id', 'candidate', 'kind', 'scaffold', 'start', 'end', 'left_anchor_start', 'right_anchor_end',
        'anchor_status', 'left_reported_placements', 'right_reported_placements', 'anchor_specificity',
        'intended_anchors_reported', 'possible_secondary_saturation',
        'hifi_status', 'hifi_source', 'hifi_source_sha256', 'primary_molecules_observed', 'molecules_bracketing_anchors',
        'molecules_passing_chain_screen', 'local_read_length_ge_anchor_span', 'longest_local_primary_read',
        'expected_bridge_count', 'support_interpretation', 'cut_authorized']
    write_table(dest/'junction_checks.tsv', fields, ({k: h.get(k, '.') for k in fields} for h in hypotheses))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--work-root', type=Path, default=Path('work'))
    p.add_argument('--hifi-results', type=Path)
    p.add_argument('--manifest', type=Path, help='JSON object: assembly ID -> ordered stage records (stage/fasta/agp/graph)')
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--minimap2', default='minimap2')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--flank', type=int, default=100000)
    p.add_argument('--anchor', type=int, default=1000)
    p.add_argument('--contact-flank', type=int, default=1000000)
    p.add_argument('--index-bases', type=int, default=8000000000)
    args = p.parse_args()
    if min(args.threads, args.flank, args.anchor, args.contact_flank, args.index_bases) <= 0 or args.anchor > args.flank:
        p.error('Require positive settings and anchor <= flank')
    args.results = args.results.resolve()
    args.out = args.out.resolve()
    args.work_root = args.work_root.resolve()
    if not args.results.is_dir():
        p.error('Results directory does not exist')
    if args.hifi_results:
        args.hifi_results = args.hifi_results.resolve()
        if not args.hifi_results.is_dir():
            p.error('HiFi results directory does not exist')
    for tool in (args.samtools, args.minimap2):
        if not shutil.which(tool):
            p.error(f'Tool not available: {tool}')
    args.out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, args.out/'trace_chimera_origins.py')
    shutil.copy2(Path(__file__).resolve().parents[2]/'py_scripts/chimera_origin.py', args.out/'chimera_origin.py')
    for tool, label in ((args.samtools, 'samtools'), (args.minimap2, 'minimap2')):
        run([tool, '--version'], args.out/(label+'.version.txt'), args.out/(label+'.version.log'))
    supplied = json.loads(args.manifest.read_text()) if args.manifest else None
    completed, scope_audit = [], []
    trace = args.results/'pipeline/pipeline_trace.txt'
    if trace.is_file():
        cached_inventory(trace, args.work_root, args.out/'cached_artifacts.tsv')
    files = sorted((args.results/'assembly/chimeras').glob('*.chimeric_joins.tsv'))
    if not files:
        raise ValueError('No called tables found')
    for path in files:
        calls = [r for r in table(path) if r.get('chromosome_member') == 'yes' and r.get('transition_lo', '.') != '.']
        audit_path = path.with_name(path.name.replace('.chimeric_joins.tsv', '.coordinate_audit.tsv'))
        audit = {r['metric']: r['value'] for r in table(audit_path)} if audit_path.is_file() else {}
        scope_audit.append(dict(table=str(path), intervals=len(calls),
            transition_assessment=audit.get('alignment_transition_assessment', 'unknown'),
            scope='existing_chromosome_candidate_intervals_only'))
        if not calls:
            continue
        assemblies = {r['assembly'] for r in calls}
        scopes = {r['join_scope'] for r in calls}
        if len(assemblies) != 1 or len(scopes) != 1 or any(r['coordinate_stage'] != 'pre_finishing' for r in calls):
            raise ValueError('Mixed assembly or coordinate stage in called table')
        assembly = next(iter(assemblies))
        if not re.fullmatch('[A-Za-z0-9_.-]+', assembly) or assembly in ('.', '..'):
            raise ValueError('Unsafe assembly ID')
        stages = supplied[assembly] if supplied is not None else discover(args.results, assembly, next(iter(scopes)))
        completed.append(process_assembly(args, assembly, calls, stages, args.out/assembly))
    write_table(args.out/'scope_audit.tsv', ['table', 'intervals', 'transition_assessment', 'scope'], scope_audit)
    report = ['# Chimera origin and junction checks', '',
              'This run performs no assembly edits. Every junction remains REVIEW.', '',
              '| Assembly | Interval | Earliest stage with an exact full-window placement |',
              '|---|---|---|']
    for result in completed:
        for key, stage in result['earliest_exact_stage'].items():
            report.append(f'| {result["assembly"]} | {key} | {stage} |')
    report.extend(['', 'Earliest verified presence is not the stage that necessarily created a join. '
        'Inspect origin_ledger.tsv for multiple placements, missing stages and Ns. '
        'A match containing Ns preserves the assembly representation, not molecular continuity.', '',
        'source_boundaries.tsv records verified coordinate projections; source_components.tsv '
        'records AGP component placement in those verified windows (input component sequence '
        'itself must be checked through the independent stage alignments). flank_placements.tsv '
        'retains exact anchor placements when the full window does not match. All base-resolved '
        'PAFs are retained, including secondary alignments. No match is not evidence of a new join.', '',
        'junction_checks.tsv screens actual N gaps, verified source boundaries, and the full '
        'uncertainty interval. Anchor specificity is a non-exhaustive search within one assessment '
        'haplotype. Read screens do not validate uniqueness or estimate expected bridge counts. '
        'Unknown physical scaffold-gap lengths remain unknown.', '',
        'alternative_partner_counts.tsv contains pooled raw Hi-C counts by scaffold. These are '
        'not normalized or independent of the data used in scaffolding. Library-specific comparisons '
        'and competitive read mapping remain follow-up tests.', '',
        'scope_audit.tsv preserves unassessed/reference states. This diagnostic consumes existing '
        'chromosome-candidate intervals; it is not an independent genome-wide chimera detector.', ''])
    (args.out/'report.md').write_text('\n'.join(report))
    (args.out/'status.json').write_text(json.dumps(dict(status='SUCCESS', assemblies=completed,
        cut_authorized=False, interpretation='Historical sequence presence and descriptive evidence only; '
        'absence of an exact match does not identify when a join arose. Anchor placements are '
        'minimap2-reported placements within one assessment haplotype, not exhaustive uniqueness. '
        'HiFi screens use retained regional records, not competitive remapping or calibrated absence tests. '
        'Hi-C partner counts are pooled, unnormalized and not independent of scaffolding.'), indent=2)+'\n')
    print(f'Review outputs: {args.out}', flush=True)


if __name__ == '__main__':
    main()
