#!/usr/bin/env python3
"""Prepare verified joins, collect descriptive evidence, and review decisions. No assembly edits."""
import argparse
from collections import defaultdict, Counter
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

from trace_chimera_origins import Fasta, run, map_queries
from chimera_origin import agp_rows, revcomp, sha, table, write_table, read_support
from junction_assessment import (validate_registry, stable_id, project_gap, choose_gaps,
                                 Lift, library_for, Contacts, decision_rows, validate_decisions)
from scan_chimera_support import scan


def write_rows(path, rows, fields=None):
    write_table(path, fields or (list(rows[0]) if rows else ['status']), rows)


def verify_agp(rows, fasta):
    ends = {}
    for r in rows:
        ends[r['object']] = r['end']
    if ends != fasta.lengths:
        raise ValueError('AGP object lengths differ from its FASTA')


def pairs_dictionary(path):
    """Read and validate the coordinate header before expensive evidence collection."""
    lengths = {}
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#chromsize:'):
                _, name, length = line.split()
                if name in lengths or int(length) <= 0:
                    raise ValueError('Duplicate or invalid pairs chromosome')
                lengths[name] = int(length)
            elif line.strip() and not line.startswith('#'):
                break
    if not lengths:
        raise ValueError('Hi-C pairs chromosome dictionary missing')
    return lengths


def find_gaps(assembly, assessment, stages, dest, samtools, scaffolds):
    """Map only identical whole scaffolds; changed stages are explicitly unresolved."""
    gaps, audits = {}, []
    current = {c: assessment.fetch(c, 0, assessment.lengths[c]) for c in scaffolds}
    signatures = defaultdict(list)
    for c, seq in current.items():
        signatures[hashlib.sha256(seq.encode()).hexdigest()].append((c, False))
        signatures[hashlib.sha256(revcomp(seq).encode()).hexdigest()].append((c, True))
    for i, spec in enumerate(stages):
        if not spec.get('agp'):
            continue
        fa, agp = Path(spec['fasta']), Path(spec['agp'])
        if not fa.is_file() or not agp.is_file():
            audits.append(dict(stage=spec['stage'], object='.', status='missing_stage_input', fasta=str(fa), agp=str(agp)))
            continue
        target = assessment if spec['stage'] == 'assessment' else Fasta(fa, dest/f'stage_{i}.fa', samtools)
        rows = agp_rows(agp)
        verify_agp(rows, target)
        by_object = defaultdict(list)
        for r in rows:
            by_object[r['object']].append(r)
        mappings = {}
        for obj, length in target.lengths.items():
            if length not in {len(s) for s in current.values()}:
                continue
            sequence = target.fetch(obj, 0, length)
            hits = [(obj, False)] if spec['stage'] == 'assessment' and obj in current else signatures.get(hashlib.sha256(sequence.encode()).hexdigest(), [])
            mappings[obj] = hits
        multiplicity = Counter(hits[0][0] for hits in mappings.values() if len(hits) == 1)
        for obj, hits in mappings.items():
            if len(hits) == 1 and multiplicity[hits[0][0]] != 1:
                hits = []
            status = 'literal_scaffold_verified' if len(hits) == 1 else 'unmapped_or_ambiguous'
            audits.append(dict(stage=spec['stage'], object=obj, status=status, fasta=str(fa), agp=str(agp)))
            if len(hits) != 1:
                continue
            scaffold, reverse = hits[0]
            length = target.lengths[obj]
            sequence = target.fetch(obj, 0, length)
            for j, r in enumerate(by_object[obj]):
                if r['kind'] not in ('N', 'U'):
                    continue
                if set(sequence[r['start']:r['end']]) != {'N'}:
                    raise ValueError('AGP gap is not an N run')
                lo, hi = project_gap(r, length, reverse)
                left = by_object[obj][j-1] if j else None
                right = by_object[obj][j+1] if j+1 < len(by_object[obj]) else None
                if reverse:
                    left, right = right, left
                key = stable_id(assembly, scaffold, lo, hi)
                if key not in gaps:
                    gaps[key] = dict(id=key, assembly=assembly, scaffold=scaffold, start=lo, end=hi,
                        origin_stage=spec['stage'], left_component=left['component'] if left else '.',
                        right_component=right['component'] if right else '.', coordinate_status=status,
                        gap_length=hi-lo, physical_gap_length='unknown')
    # Per-scaffold status makes missing/changed earlier sequences visible even if no lengths match.
    for c in scaffolds:
        audits.append(dict(stage='origin_summary', object=c,
                           status='Earlier stages require exact whole-scaffold correspondence; absent recovery is not evidence of a round2 origin.', fasta='.', agp='.'))
    write_rows(dest/'stage_projection_audit.tsv', audits)
    return list(gaps.values())


def bam_frame(meta, assembly, digest, bam, assessment, samtools):
    if meta.get('assembly') != assembly or meta.get('sha256') != digest or meta.get('coordinate_stage') != 'pre_finishing' or not meta.get('hifi'):
        raise ValueError('HiFi provenance disagrees with assessment')
    run([samtools, 'quickcheck', bam])
    header = subprocess.check_output([samtools, 'view', '-H', str(bam)], text=True)
    lengths = {}
    for line in header.splitlines():
        if line.startswith('@SQ\t'):
            tags = dict(f.split(':', 1) for f in line.split('\t')[1:])
            if tags['SN'] in lengths:
                raise ValueError('Duplicate BAM sequence')
            lengths[tags['SN']] = int(tags['LN'])
    if lengths != assessment.lengths:
        raise ValueError('HiFi BAM coordinate dictionary mismatch')


def pipe_fastq(samtools, bam, names, fastq, log):
    with fastq.open('w') as out, log.open('w') as err:
        upstream = subprocess.Popen([samtools, 'view', '-b', '-N', str(names), '-F', '2304', str(bam)], stdout=subprocess.PIPE, stderr=err)
        try:
            downstream = subprocess.run([samtools, 'fastq', '-n', '-'], stdin=upstream.stdout, stdout=out, stderr=err)
        finally:
            upstream.stdout.close()
        rc = upstream.wait()
        if rc or downstream.returncode:
            raise RuntimeError('Read extraction failed')


def hifi_evidence(a, job, joint, names_by_scaffold, dest):
    assembly, assessment, digest = job['assembly'], job['fasta'], job['digest']
    context = a.hifi_results/'assembly/chimeras/sequence_context'/f'{assembly}.sequence_context'
    bam, provenance = context/'hifi.bam', context/'provenance.json'
    if not bam.is_file() or not provenance.is_file():
        return dict(status='unavailable', reason='retained full HiFi BAM or provenance missing')
    meta = json.loads(provenance.read_text())
    bam_frame(meta, assembly, digest, bam, assessment, a.samtools)
    shutil.copy2(provenance, dest/'hifi_source_provenance.json')
    regions = [dict(id=s['id'], scaffold=s['scaffold'], start=int(s['start']), end=int(s['end']), role=s['kind']) for s in job['seeds']]
    regions += job['selected']
    reads = set()
    for r in regions:
        lo, hi = max(0, int(r['start'])-a.flank), min(assessment.lengths[r['scaffold']], int(r['end'])+a.flank)
        sam = dest/f"{r['id']}.original.sam"
        run([a.samtools, 'view', '-h', bam, f"{r['scaffold']}:{lo+1}-{hi}"], stdout=sam)
        with sam.open() as handle:
            reads.update(line.split('\t', 1)[0] for line in handle if not line.startswith('@'))
    names = dest/'selected_read_names.txt'
    names.write_text(''.join(n+'\n' for n in sorted(reads)))
    fastq, sam = dest/'selected.fastq', dest/'competitive.sam'
    if reads:
        pipe_fastq(a.samtools, bam, names, fastq, dest/'read_extraction.log')
        if fastq.stat().st_size == 0:
            raise ValueError('Selected reads yielded no primary FASTQ sequences')
        cmd = [a.minimap2, '-ax', 'map-hifi', '-I', str(a.index_bases), '--secondary=yes', '-N', '50',
               '-t', str(a.threads), joint, fastq]
        (dest/'competitive.command.json').write_text(json.dumps([str(x) for x in cmd]))
        run(cmd, stdout=sam, stderr=dest/'competitive.log')
    else:
        sam.write_text('')
    # Preserve all emitted placements without SEQ/QUAL in the transferable packet.
    with sam.open() as src, gzip.open(dest/'competitive_alignments.tsv.gz', 'wt') as out:
        out.write('qname\tflag\trname\tpos_1based\tmapq\tcigar\ttags\n')
        for line in src:
            if line.startswith('@'):
                continue
            f = line.rstrip().split('\t')
            out.write('\t'.join(f[:6]+[';'.join(f[11:])])+'\n')
    support = []
    for r in regions:
        scaffold, start, end = r['scaffold'], int(r['start']), int(r['end'])
        lo, hi = max(0, start-a.flank), min(assessment.lengths[scaffold], end+a.flank)
        renamed = names_by_scaffold[assembly, scaffold]
        local = dest/f"{r['id']}.competitive.sam"
        # Keep only overlapping records; avoid counting distant same-scaffold reads as local opportunity.
        with sam.open() as src, local.open('w') as out:
            from scan_chimera_support import alignment
            for line in src:
                if line.startswith('@'):
                    continue
                f = line.rstrip().split('\t')
                if f[2] != renamed or f[5] == '*':
                    continue
                s, e, _, _ = alignment(f)
                if s < hi and lo < e:
                    out.write(line)
        sequence = assessment.fetch(scaffold, lo, hi)
        _, track = scan(local, dict(scaffold=renamed, start=lo, end=hi, lo=start, hi=end), sequence, step=a.step, anchor=a.anchor)
        write_rows(dest/f"{r['id']}.support_track.tsv", track)
        astart, aend = start-a.anchor, end+a.anchor
        valid = astart >= 0 and aend <= assessment.lengths[scaffold]
        if valid:
            valid = set(assessment.fetch(scaffold, astart, start)+assessment.fetch(scaffold, end, aend)) <= set('ACGT')
        row = dict(id=r['id'], role=r['role'], scaffold=scaffold, start=start, end=end,
                   anchor_status='sequence_available' if valid else 'terminal_or_ambiguous_bases',
                   selection_bias='reads selected from existing same-haplotype mapping; unmapped reads are not assessed')
        if valid:
            with local.open() as handle:
                row.update(read_support(handle, renamed, astart, aend))
        support.append(row)
    fields = list(dict.fromkeys(k for row in support for k in row))
    write_rows(dest/'hifi_support.tsv', support, fields)
    return dict(status='measured', selected_read_names=len(reads), bam=str(bam.resolve()),
                bam_dictionary_verified=True, limitation='Competitive against supplied current assemblies only. Shared haplotypes may lower MAPQ. No calibrated absence or phase-consistency verdict.')


def comparison_context(a, jobs):
    """Base-resolved context for every requested interval and candidate/control gap."""
    dest = a.out/'sequence_comparisons'
    dest.mkdir()
    queries, records = dest/'windows.fa', []
    with queries.open('w') as out:
        for job in jobs:
            for r in job['seeds']+job['selected']:
                lo=max(0,int(r['start'])-a.flank)
                hi=min(job['fasta'].lengths[r['scaffold']],int(r['end'])+a.flank)
                out.write(f">{r['id']}\n{job['fasta'].fetch(r['scaffold'],lo,hi)}\n")
                records.append(dict(id=r['id'],assembly=job['assembly'],scaffold=r['scaffold'],start=lo,end=hi,
                                    core_start=int(r['start']),core_end=int(r['end'])))
    write_rows(dest/'windows.tsv',records)
    paths = []
    if a.peer_fasta_dir:
        paths += [('current_peer',str(p.resolve())) for p in sorted(a.peer_fasta_dir.iterdir())
                  if p.is_file() and p.suffix in ('.fa','.fasta','.fna')]
    if a.history:
        meta=json.loads((a.history/'provenance.json').read_text())
        for r in meta['inputs']:
            if r['role']=='historical':
                path=Path(r['path'])
                if sha(path)!=r['sha256']:
                    raise ValueError('Historical FASTA has changed since comparison')
                paths.append((r['id'],str(path.resolve())))
    audit=[]
    seen=set()
    for i,(label,filename) in enumerate(paths):
        if filename in seen:
            continue
        seen.add(filename)
        print(f'Comparing junction windows against {filename}',flush=True)
        target=Fasta(Path(filename),dest/f'peer_{i}.fa',a.samtools)
        hits=map_queries(a,target,queries,dest/f'peer_{i}_windows')
        summary=[]
        for h in hits:
            summary.append(dict(query=h['query'],target=h['target'],query_start=h['query_start'],query_end=h['query_end'],
                                target_start=h['target_start'],target_end=h['target_end'],strand=h['strand'],mapq=h['mapq'],
                                matches=h['matches'],alignment_columns=h['block']))
        write_rows(dest/f'peer_{i}_alignments.tsv',summary)
        audit.append(dict(label=label,path=filename,sha256=sha(Path(filename)),alignment_file=f'peer_{i}_windows.paf'))
    write_rows(dest/'sources.tsv',audit)
    return len(audit)


def graph_context(a, jobs):
    """Retain native read-placement/topology records for implicated raw contigs."""
    for job in jobs:
        ledger=table(a.origins/job['assembly']/'origin_ledger.tsv')
        targets={r['target'] for r in ledger if r['stage']=='raw_contigs' and r['target']!='.'}
        targets.update(c for g in job['selected'] for c in (g['left_component'],g['right_component']) if c!='.')
        graph=next((s.get('graph') for s in job['stages'] if s['stage']=='raw_contigs'),None)
        if not graph or not Path(graph).is_file():
            (job['dest']/'graph_status.json').write_text(json.dumps(dict(status='unavailable')))
            continue
        counts=Counter()
        with Path(graph).open() as src, gzip.open(job['dest']/'raw_graph_context.gfa.gz','wt') as out:
            for line in src:
                f=line.rstrip().split('\t')
                keep=(f[0] in ('A','S') and len(f)>1 and f[1] in targets) or (f[0]=='L' and len(f)>3 and (f[1] in targets or f[3] in targets))
                if not keep:
                    continue
                counts[f[0]]+=1
                if f[0]=='S' and len(f)>2:
                    f[2]='*'
                out.write('\t'.join(f)+'\n')
        (job['dest']/'graph_status.json').write_text(json.dumps(dict(status='native_records_collected',path=graph,
            targets=sorted(targets),record_counts=dict(counts),interpretation='Native A/L records, no inferred read path, overlap or biological validation. Purged-name aliases are not guessed.'),indent=2))


def hic_evidence(a, job, dest):
    assembly, specs, assessment = job['assembly'], job['stages'], job['fasta']
    final = next(s for s in specs if s['stage'] == 'assessment')
    is_round2 = any(s['stage'] == 'scaffold_corrected' for s in specs)
    stage = 'scaffold' if is_round2 else 'contig'
    base = a.results/'bam/hic'/stage
    pairs, manifest = base/'filtered'/f'{assembly}.pairs.gz', base/'raw'/f'{assembly}.readsets.tsv'
    if not pairs.is_file() or not manifest.is_file():
        return dict(status='unavailable', reason='filtered source pairs or read-set manifest missing', pairs=str(pairs), manifest=str(manifest))
    rows = table(manifest)
    prefixes = [r['qname_prefix'] for r in rows]
    if not rows or any(not p for p in prefixes) or any(p.startswith(q) for i, p in enumerate(prefixes) for j, q in enumerate(prefixes) if i != j):
        raise ValueError('Missing, duplicate or overlapping read-set prefixes')
    agp = agp_rows(final['agp'])
    verify_agp(agp, assessment)
    lift = Lift(agp)
    source_stage = 'scaffold_corrected' if is_round2 else 'round1_input'
    override = getattr(a, 'hic_input_map', {}).get(assembly)
    source = next((s for s in specs if s['stage'] == source_stage), None)
    if not override and (not source or not Path(source['fasta']).is_file()):
        raise ValueError('Cannot identify exact Hi-C mapping input; supply --hic-inputs')
    source_path = Path(override['fasta']) if override else Path(source['fasta'])
    if override and sha(source_path) != override['sha256']:
        raise ValueError('Explicit Hi-C input FASTA checksum mismatch')
    sourcefa = Fasta(source_path, dest/'hic_input.fa', a.samtools)
    if pairs_dictionary(pairs) != sourcefa.lengths:
        raise ValueError('Hi-C pairs dictionary differs from input FASTA; supply --hic-inputs with the exact last-scaffolding mapping input')
    for r in agp:
        if r['kind'] not in ('N', 'U') and (r['component'] not in sourcefa.lengths or r['component_end'] > sourcefa.lengths[r['component']]):
            raise ValueError('Hi-C source FASTA does not match AGP components')
    contacts = Contacts(job['selected'], assessment.lengths, a.contact_flank)
    audits = Counter()
    with gzip.open(pairs, 'rt') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip().split('\t')
            if len(f) < 8:
                raise ValueError('Malformed pairs')
            audits['input_pairs'] += 1
            if f[7] != 'UU':
                continue
            lib = library_for(f[0], rows)
            audits['UU_'+lib] += 1
            x, y = lift.locate(f[1], int(f[2])-1), lift.locate(f[3], int(f[4])-1)
            if x is None or y is None:
                audits['unmapped_or_ambiguous_pairs'] += 1
                continue
            contacts.add(lib, x, y)
    output, alternatives = [], []
    libraries = sorted({r['library_id'] for r in rows} | set(contacts.total))
    for lib in libraries:
        for g in job['selected']:
            output.append(dict(id=g['id'], role=g['role'], library=lib, crossing_pairs=contacts.cross[lib, g['id']],
                projected_UU_pairs=contacts.total[lib], crossing_per_million_projected_UU=contacts.cross[lib, g['id']]*1e6/contacts.total[lib] if contacts.total[lib] else '.',
                left_window_bp=contacts.windows[g['id'], 'left'][2]-contacts.windows[g['id'], 'left'][1],
                right_window_bp=contacts.windows[g['id'], 'right'][2]-contacts.windows[g['id'], 'right'][1],
                left_external_contacts=contacts.margin[lib,g['id'],'left'], right_external_contacts=contacts.margin[lib,g['id'],'right']))
            for side in ('left', 'right'):
                partners = [(k, v) for k, v in contacts.partners.items() if k[:3] == (lib,g['id'],side)]
                for k, count in sorted(partners, key=lambda kv: (-kv[1], kv[0]))[:10]:
                    c, b = k[3:]
                    alternatives.append(dict(id=g['id'], library=lib, side=side, partner_scaffold=c,
                        partner_start=b*a.contact_flank, partner_end=min((b+1)*a.contact_flank, assessment.lengths[c]),
                        contacts=count, fraction_of_flank_external_contacts=count/contacts.margin[lib,g['id'],side]))
    write_rows(dest/'hic_by_library.tsv', output)
    write_rows(dest/'hic_alternative_bins.tsv', alternatives)
    shutil.copy2(manifest, dest/'hic_readsets.tsv')
    return dict(status='measured', audits=dict(audits), pairs=str(pairs.resolve()),
                readsets_sha256=sha(manifest), source_fasta=str(source_path.resolve()),
                source_fasta_sha256=sha(source_path), agp_sha256=sha(Path(final['agp'])),
                coordinate_basis='Published pipeline layout plus complete source dictionary and AGP. No sequence checksum embedded in pairs.',
                limitation='Same scaffolding reads; descriptive counts, not independent validation or a distance/mappability-adjusted score.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registry', type=Path, required=True)
    p.add_argument('--origins', type=Path, required=True)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--hifi-results', type=Path)
    p.add_argument('--hic-inputs', type=Path, help='Optional TSV: assembly, fasta, sha256; exact FASTAs used for last-scaffolding Hi-C mapping')
    p.add_argument('--history', type=Path, help='Historical comparison directory with provenance.json')
    p.add_argument('--peer-fasta-dir', type=Path, help='Explicit directory of comparable current pre-finishing FASTAs')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--prepare-only', action='store_true')
    p.add_argument('--decisions', type=Path, help='Validate a reviewed table against this new coordinate-verified packet; no cuts exported')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--flank', type=int, default=100000)
    p.add_argument('--anchor', type=int, default=1000)
    p.add_argument('--step', type=int, default=1000)
    p.add_argument('--contact-flank', type=int, default=100000)
    p.add_argument('--control-count', type=int, default=3)
    p.add_argument('--control-radius', type=int, default=1000000)
    p.add_argument('--index-bases', type=int, default=8000000000)
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--minimap2', default='minimap2')
    a = p.parse_args()
    if min(a.threads,a.flank,a.anchor,a.step,a.contact_flank,a.control_radius,a.index_bases) <= 0 or a.control_count < 0 or a.anchor > a.flank:
        p.error('Invalid resource/window settings')
    seeds = validate_registry(table(a.registry))
    if not seeds:
        p.error('Empty registry')
    if len({s['sample'] for s in seeds}) != 1:
        p.error('Use one individual per competitive-mapping packet')
    a.out.mkdir(parents=True, exist_ok=False)
    (a.out/'parameters.json').write_text(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},indent=2))
    a.hic_input_map = {}
    if a.hic_inputs:
        for r in table(a.hic_inputs):
            if r['assembly'] in a.hic_input_map or r['assembly'] not in {s['assembly'] for s in seeds}:
                raise ValueError('Unknown or duplicate assembly in Hi-C input manifest')
            a.hic_input_map[r['assembly']] = r
        shutil.copy2(a.hic_inputs, a.out/'hic_input_manifest.tsv')
    shutil.copy2(a.registry, a.out/'registry.tsv')
    for tool, flag in ((a.samtools,'--version'),(a.minimap2,'--version')):
        run([tool,flag], stdout=a.out/(Path(tool).name+'.version.txt'))
    for src in [Path(__file__), Path(__file__).with_name('trace_chimera_origins.py'), Path(__file__).with_name('scan_chimera_support.py'),
                Path(__file__).parents[2]/'py_scripts/chimera_origin.py', Path(__file__).parents[2]/'py_scripts/junction_assessment.py']:
        shutil.copy2(src,a.out/src.name)
    jobs, allgaps, links, provenance = [], [], [], []
    for i, assembly in enumerate(sorted({s['assembly'] for s in seeds})):
        dest = a.out/f'assembly_{i}'
        dest.mkdir()
        selected_seeds = [s for s in seeds if s['assembly'] == assembly]
        origin = json.loads((a.origins/assembly/'provenance.json').read_text())
        shutil.copy2(a.origins/assembly/'provenance.json',dest/'origin_provenance.json')
        stages = origin['stages']
        spec = next(s for s in stages if s['stage'] == 'assessment')
        write_rows(dest/'stage_checksums.tsv',[dict(stage=s['stage'],fasta=s['fasta'],
            agp=s.get('agp',''),agp_sha256=sha(Path(s['agp'])) if s.get('agp') and Path(s['agp']).is_file() else '.') for s in stages])
        fa = Path(spec['fasta'])
        digest = sha(fa)
        if {s['assessment_sha256'] for s in selected_seeds} != {digest} or origin['assessment_sha256'] != digest:
            raise ValueError('Registry/origin/assessment checksum mismatch')
        fasta = Fasta(fa,dest/'assessment.fa',a.samtools)
        for s in selected_seeds:
            if s['scaffold'] not in fasta.lengths or int(s['end']) > fasta.lengths[s['scaffold']]:
                raise ValueError('Registry interval outside assessment')
        gaps = find_gaps(assembly,fasta,stages,dest,a.samtools,{s['scaffold'] for s in selected_seeds})
        selected, mapping = choose_gaps(selected_seeds,gaps,a.control_count,a.control_radius)
        # Only generic physical-gap IDs appear in algorithms; investigation IDs are data.
        write_rows(dest/'all_verified_gaps.tsv',gaps)
        write_rows(dest/'selected_gaps.tsv',selected)
        allgaps += gaps
        links += mapping
        jobs.append(dict(assembly=assembly,fasta=fasta,digest=digest,stages=stages,seeds=selected_seeds,selected=selected,dest=dest))
    write_rows(a.out/'seed_gap_links.tsv',links)
    write_rows(a.out/'decisions.tsv',decision_rows(seeds,links))
    graph_context(a,jobs)
    if a.decisions:
        validate_decisions(table(a.decisions),seeds,allgaps)
        shutil.copy2(a.decisions,a.out/'reviewed_decisions.tsv')
    # Prefix references so haplotypes compete without ambiguous scaffold identifiers.
    joint, names = a.out/'competitive_reference.fa', {}
    if sum(sum(j['fasta'].lengths.values()) for j in jobs) >= a.index_bases:
        raise ValueError('Combined reference exceeds single-index budget')
    if not a.prepare_only and a.hifi_results:
        # Never competitively combine unrelated individuals without an explicit study design.
        samples = set()
        for j in jobs:
            meta_path=a.hifi_results/'assembly/chimeras/sequence_context'/f"{j['assembly']}.sequence_context/provenance.json"
            if meta_path.is_file():
                samples.add(json.loads(meta_path.read_text())['sample'])
        if samples and samples != {s['sample'] for s in seeds}:
            raise ValueError('Run a separate packet per individual for competitive HiFi mapping')
        with joint.open('w') as out:
            for i, j in enumerate(jobs):
                with j['fasta'].path.open() as source:
                    for line in source:
                        if line.startswith('>'):
                            name=line[1:].split()[0]
                            names[j['assembly'],name]=f'a{i}__{name}'
                            out.write('>'+names[j['assembly'],name]+'\n')
                        else:
                            out.write(line)
        write_rows(a.out/'competitive_reference_names.tsv',[dict(assembly=k[0],scaffold=k[1],mapped_name=v) for k,v in names.items()])
    for j in jobs:
        print(f"Assessing {j['assembly']}: {len(j['selected'])} candidate/control gaps",flush=True)
        hic = hic_evidence(a,j,j['dest']) if not a.prepare_only else dict(status='not_requested')
        hifi = hifi_evidence(a,j,joint,names,j['dest']) if not a.prepare_only and a.hifi_results else dict(status='not_requested')
        provenance.append(dict(assembly=j['assembly'],assessment_sha256=j['digest'],folder=j['dest'].name,hifi=hifi,hic=hic))
    comparisons=comparison_context(a,jobs) if not a.prepare_only and (a.history or a.peer_fasta_dir) else 0
    (a.out/'status.json').write_text(json.dumps(dict(status='SUCCESS',cut_authorized=False,assemblies=provenance,comparison_assemblies=comparisons),indent=2))
    (a.out/'report.md').write_text('# Junction assessment packet\n\nNo assembly edits or cutting instructions are produced. See status.json for evidence availability. decisions.tsv is a review worksheet, not an automatic classifier. Nearby controls are unvalidated comparison joins. Gap coordinates are literal verified AGP gaps; changed earlier scaffolds remain unresolved. Competitive mapping uses only selected reads originally mapped near the requested regions and the supplied current assemblies. No calibrated absence test, independent Hi-C validation or biological-fusion verdict is made.\n')


if __name__ == '__main__':
    main()
