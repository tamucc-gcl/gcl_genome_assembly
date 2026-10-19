#!/usr/bin/env python3
"""Read-only junction investigation; inputs and hypotheses are supplied in a manifest.

This deliberately does not change production discovery, scoring, or assemblies.
All intervals are zero-based, half-open; a probe is the number of bases before a cut.
"""
import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import statistics
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'py_scripts'))
from chimera_controls import agp, literal_gaps, molecule_chains, qualified_molecule
from chimera_tracks import blocks, LabelIndex, position
from misassembly_core import read_paf, independent_labels, sha

COMMANDS=[]


def table(path):
    with open(path, encoding='utf-8') as src:
        return list(csv.DictReader(src, delimiter='\t'))


def write(path, rows, fields):
    with open(path, 'w', newline='', encoding='utf-8') as dst:
        w = csv.DictWriter(dst, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def save(path, data):
    path=Path(path);text=json.dumps(data, indent=2)
    if not path.exists() or path.read_text(encoding='utf-8')!=text:
        path.write_text(text, encoding='utf-8')


def run(command, output=None):
    COMMANDS.append(list(map(str,command)))
    print('[join investigation] ' + ' '.join(map(str, command[:4])), file=sys.stderr, flush=True)
    if output:
        with open(output, 'w') as dst:
            subprocess.run(list(map(str, command)), stdout=dst, check=True)
    else:
        subprocess.run(list(map(str, command)), check=True)


def command_args(task, executable):
    """Read the actual cached command, never infer an input from a similar filename."""
    tokens = shlex.split((task / '.command.sh').read_text())
    index = next((i for i, t in enumerate(tokens) if Path(t).name == executable), None)
    if index is None:
        raise ValueError(f'{executable} absent from {task}/.command.sh')
    result = {}
    for i in range(index + 1, len(tokens)-1):
        if tokens[i].startswith('--'):
            result[tokens[i][2:]] = tokens[i+1]
    return result


def required(path):
    path = Path(path).resolve()
    if not path.is_file():
        raise ValueError(f'Required cached input missing: {path}')
    return str(path)


def stage_commands(task):
    """Read only nxf_stage, excluding similarly named output-copy commands."""
    wrapper=task/'.command.run'
    if not wrapper.is_file():return []
    commands=[];inside=False
    for line in wrapper.read_text().splitlines():
        if re.match(r'^\s*nxf_stage\(\)\s*\{\s*$',line):inside=True;continue
        if inside and re.match(r'^\s*}\s*$',line):break
        if inside and re.match(r'^\s*(?:cp|ln)\s',line):commands.append(line.strip())
    return commands


def staged_links(task):
    """Read literal input source/destination mappings for copy and link modes."""
    result=[]
    for line in stage_commands(task):
        try:tokens=shlex.split(line,comments=True)
        except ValueError:continue
        operands=[t for t in tokens[1:] if not t.startswith('-')]
        if len(operands)!=2:continue
        source,destination=operands
        # Nextflow writes concrete paths. Expanding shell expressions here would
        # turn metadata parsing into shell execution or an unsupported guess.
        if any(c in source+destination for c in ('$','`','\n')):continue
        result.append((Path(os.path.abspath(task/destination)),Path(os.path.abspath(task/source))))
    return result


class RetainedInputError(ValueError):
    def __init__(self,requested,attempts,wrappers):
        self.details=dict(requested=str(requested),attempted_paths=attempts,wrappers=wrappers)
        super().__init__(f'Required retained input unavailable: {requested}. '
                         'No existing source was recoverable from recorded Nextflow input staging.')


def retained_input(task,value,hint=None):
    """Recover a cleaned scratch-stage input through its exact recorded link chain.

    A staged directory link can supply a file beneath it. Upstream staging links
    can also have been removed; follow their wrappers recursively. No work-tree
    search by basename, copied input, or recreated symlink is used.
    """
    requested=Path(os.path.abspath(task/value));seen=set();attempts=[];wrappers={}
    def resolve(path):
        path=Path(os.path.abspath(path))
        if path.is_file():return str(path.resolve())
        if path in seen or len(seen)>=64:return None
        seen.add(path);attempts.append(str(path))
        for parent in path.parents:
            matches=[]
            if (parent/'.command.run').is_file() and str(parent) not in wrappers:
                wrappers[str(parent)]=dict(path=str(parent/'.command.run'),commands=stage_commands(parent),
                    mappings=[dict(destination=str(d),source=str(s),source_exists=s.is_file() or s.is_dir())
                              for d,s in staged_links(parent)])
            for destination,source in staged_links(parent):
                try:rest=path.relative_to(destination)
                except ValueError:continue
                matches.append((len(destination.parts),source/rest))
            for _,source in sorted(matches,key=lambda item:item[0],reverse=True):
                found=resolve(source)
                if found:return found
        return None
    # Absolute retained paths already recorded in the evidence take precedence
    # over transient staging names; downstream checksum checks still apply.
    if hint:
        found=resolve(Path(hint))
        if found:return found
    found=resolve(requested)
    if found:return found
    raise RetainedInputError(requested,attempts,list(wrappers.values()))


def cached_task(evidence):
    # mapping_manifest provenance is a cached <assembly>.evidence/evidence.json.
    return Path(evidence).parent.parent


def windows(event, size, padding=2_000_000):
    sites = [dict(id=event['id']+'T', kind='transition_midpoint', lo=event['lo'],
                  hi=event['hi'], point=(event['lo']+event['hi'])//2)]
    for source, kind in ((event.get('alternatives', []), 'regional_trough'),
                         (event.get('cut_options', []), 'gap_option')):
        for s in source:
            point = s.get('probe_bp') if kind == 'regional_trough' else s.get('cut_bp')
            if point is None: point = (s['lo']+s['hi'])//2
            sites.append(dict(id=s['id'], kind=kind, lo=s['lo'], hi=s['hi'], point=point))
    lo = max(0, min(s['lo'] for s in sites)-padding)
    hi = min(size, max(s['hi'] for s in sites)+padding)
    return sites, lo, hi


def prepare(a):
    """Resolve exact retained task inputs from the published packet and its cache."""
    a.out.mkdir(parents=True, exist_ok=True)
    source = a.source.resolve()
    mappings = json.loads((source/'mapping_manifest.json').read_text())
    selected = table(a.focus)
    if len({(r['assembly'], r['event']) for r in selected}) != len(selected):
        raise ValueError('Duplicate investigation target')
    grouped = defaultdict(list)
    for row in selected: grouped[row['assembly']].append(row)
    jobs = [];failures=[];resolution=dict(assemblies=[])
    def resolve_input(audit,key,task,value,hint=None):
        try:
            path=retained_input(task,value,hint)
            audit['inputs'][key]=dict(status='resolved',path=path)
            return path
        except RetainedInputError as error:
            audit['inputs'][key]=dict(status='unavailable',**error.details)
            failures.append(audit['assembly']+': '+key)
            return None
        finally:save(a.out/'input_resolution.json',resolution)
    for assembly, rows in grouped.items():
        published = source/'assemblies'/assembly/'evidence.json'
        data = json.loads(published.read_text())
        task = cached_task(mappings[assembly]['provenance'])
        args = command_args(task, 'misassembly_assess.py')
        audit=dict(assembly=assembly,assessment_task=str(task),inputs={});resolution['assemblies'].append(audit)
        hints={'fasta':data.get('assessment_fasta') or mappings[assembly].get('fasta')}
        paths = {key: resolve_input(audit,key,task,value,hints.get(key)) for key, value in args.items()
                 if key in ('fasta', 'discovery', 'mapping', 'agp', 'pairs', 'source', 'libraries')}
        bam=resolve_input(audit,'bam',task,mappings[assembly]['bam'])
        index=resolve_input(audit,'index',task,mappings[assembly]['index'])
        for key in ('fasta', 'discovery', 'mapping', 'agp', 'pairs', 'source', 'libraries'):
            if key not in paths: raise ValueError(f'{assembly}: {key} unavailable in retained task')
        for key, expected in (('fasta', data['assessment_sha256']),
                              ('agp', data['hic_audit']['agp_sha256']),
                              ('source', data['hic_audit']['source_sha256'])):
            if paths[key] is not None and sha(paths[key]) != expected: raise ValueError(f'{assembly}: changed {key}')
        if data.get('coordinate_stage') != 'pre_finishing':
            raise ValueError('Expected explicit pre_finishing evidence coordinates')
        if paths['discovery'] is None:continue
        # The staged discovery path resolves back to its actual Nextflow task.
        discovery_task = Path(paths['discovery']).parent.parent
        discovery_args = command_args(discovery_task, 'misassembly_discover.py')
        catalog = resolve_input(audit,'catalog',discovery_task,discovery_args['catalog'])
        alignments = json.loads(discovery_args['alignments'])
        for pair in alignments:
            pair['path'] = resolve_input(audit,'alignment:'+pair['a']+':'+pair['b'],discovery_task,pair['path'])
        if (catalog is None or bam is None or index is None or
                any(path is None for path in paths.values()) or any(pair['path'] is None for pair in alignments)):continue
        if { (r['a'], r['b'], r['sha256']) for r in data['alignments'] } != {
                (r['a'], r['b'], sha(r['path'])) for r in alignments }:
            raise ValueError('Cohort alignment provenance differs from reviewed evidence')
        mapping = json.loads(Path(paths['mapping']).read_text())
        if (mapping['assessment_sha256'] != data['assessment_sha256'] or
                mapping.get('assembly')!=assembly or mapping.get('sample')!=data['sample'] or
                mapping.get('coordinate_stage')!='pre_finishing'):
            raise ValueError('BAM provenance differs from assessment sequence')
        if 'mapping' in data and Path(mappings[assembly]['bam']).resolve()!=Path(data['mapping']['bam']).resolve():
            raise ValueError('Published BAM manifest differs from evidence provenance')
        run(['samtools', 'quickcheck', bam])
        sizes = {f[0]:int(f[1]) for f in (line.split('\t') for line in
                 subprocess.check_output(['samtools','idxstats',bam],text=True).splitlines()) if f[0]!='*'}
        fasta_sizes = {}
        # Do not write indexes beside immutable cache inputs.
        with open(paths['fasta']) as handle:
            name = None
            for line in handle:
                if line.startswith('>'): name=line[1:].split()[0]; fasta_sizes[name]=0
                else: fasta_sizes[name] += len(line.strip())
        if sizes != fasta_sizes: raise ValueError('BAM/FASTA sequence dictionaries differ')
        events=[]
        for row in rows:
            event=next(e for e in data['events'] if e['id']==row['event'])
            sites,lo,hi=windows(event,sizes[event['scaffold']])
            events.append(dict(event=event['id'],role=row['role'],scaffold=event['scaffold'],
                               left=event['left'],right=event['right'],lo=lo,hi=hi,sites=sites,
                               transition_lo=event['lo'],transition_hi=event['hi']))
        final_map = required(a.final_dir/(assembly+'.name_map.tsv'))
        # Names after finishing are context, not an exact coordinate lift: end trimming
        # and extension can change positions. Never silently add a parent offset.
        jobs.append(dict(assembly=assembly,sample=data['sample'],sizes=sizes,events=events,
                         evidence=required(published),catalog=catalog,alignments=alignments,
                         bam=bam,index=index,final_map=final_map,paths=paths,
                         source_evidence_sha256=sha(published),catalog_sha256=sha(catalog),final_map_sha256=sha(final_map),
                         input_stats={str(p):dict(size=Path(p).stat().st_size,mtime_ns=Path(p).stat().st_mtime_ns)
                                      for p in [bam,index,*paths.values()]},
                         assessment_sha256=data['assessment_sha256'],
                         coordinate_stage='pre_finishing',
                         hic_expected=data['hic_audit']['audit']['input_pairs']))
    if failures:
        raise ValueError('Unresolved investigation inputs: '+', '.join(failures)+'. '
                         f'Details: {a.out / "input_resolution.json"} (included in the review archive).')
    if len({j['sample'] for j in jobs}) != 1:
        raise ValueError('This competitive read investigation expects haplotypes of one individual')
    for job in jobs:
        job['competitive_fastas']={j['assembly']:j['paths']['fasta'] for j in jobs}
        save(a.out/(job['assembly']+'.json'),job)
    save(a.out/'manifest.json',dict(jobs=jobs,source=str(source),focus=selected,
                                  coordinate_convention='zero-based half-open; cut_bp bases before cut'))


def local_hits(pair, job, regions):
    result=defaultdict(list)
    for f in read_paf(pair['path'],job['assembly']==pair['b'],set(regions)):
        if any(int(f[2])<hi and int(f[3])>lo for lo,hi in regions[f[0]]): result[f[0]].append(f)
    return result


def union(intervals):
    result=[]
    for lo,hi in sorted(intervals):
        if hi<=lo: continue
        if result and lo<=result[-1][1]: result[-1]=(result[-1][0],max(hi,result[-1][1]))
        else: result.append((lo,hi))
    return result


def overlap(intervals,lo,hi):
    return sum(max(0,min(hi,b)-max(lo,a)) for a,b in intervals)


def consensus_masks(observations, focal_sample, minimum=2):
    """Unique-locus support, two independent individuals, no informative disagreement."""
    events=defaultdict(list)
    for sample,lo,hi,chrom in observations:
        if sample==focal_sample or not chrom: continue
        events[lo].append((sample,chrom,1));events[hi].append((sample,chrom,-1))
    active=defaultdict(Counter);points=sorted(events);out=[]
    for i,lo in enumerate(points[:-1]):
        for sample,chrom,change in events[lo]:
            active[sample][chrom]+=change
            if not active[sample][chrom]:del active[sample][chrom]
        votes=[next(iter(c)) for c in active.values() if len(c)==1]
        conflict=any(len(c)>1 for c in active.values())
        if not conflict and len(votes)>=minimum and len(set(votes))==1:
            chrom=votes[0];hi=points[i+1]
            if out and out[-1]['hi']==lo and out[-1]['chrom']==chrom:out[-1]['hi']=hi
            else:out.append(dict(lo=lo,hi=hi,chrom=chrom))
    return out


def local_tracks(hits, labels, lo, hi, mode):
    """Clip CIGAR blocks before the interval sweep, including reverse-strand offsets."""
    index=LabelIndex(labels);segments=[];edges=defaultdict(list)
    for f in hits:
        for b in blocks(f,compact=True):
            start=max(lo,b['lo']);end=min(hi,b['hi'])
            if end<=start:continue
            copy=dict(b,lo=start,hi=end,
                      t=b['t']+(start-b['lo'] if b['strand']=='+' else b['hi']-end))
            i=len(segments);segments.append(copy);edges[start].append((i,1));edges[end].append((i,-1))
            for edge in index.inside(copy['target'],copy['t'],copy['t']+end-start,True):
                edges[start+edge-copy['t'] if copy['strand']=='+' else end-(edge-copy['t'])]
    active=set();points=sorted(edges);out=[]
    for j,start in enumerate(points[:-1]):
        for i,delta in edges[start]:
            if delta==1:active.add(i)
            else:active.discard(i)
        end=points[j+1];chrom=None;status='unaligned_or_filtered'
        good=[segments[i] for i in active if segments[i]['mapq']>=20 and segments[i]['identity']>=.9 and segments[i]['target'] in labels]
        if good:
            best=max(good,key=lambda b:(b['mapq'],b['identity']));mid=(start+end)//2
            chrom=index.at(best['target'],position(best,mid))
            competing=[segments[i] for i in active if segments[i]['identity']>=.95*best['identity'] and
                       (segments[i]['target']!=best['target'] or segments[i]['strand']!=best['strand'] or abs(position(segments[i],mid)-position(best,mid))>1000)]
            if competing and (mode=='locus_unique' or any(index.at(b['target'],position(b,mid))!=chrom for b in competing)):
                chrom=None;status='ambiguous'
            else:status=('chromosome_assigned_locus_ambiguous' if competing else 'assigned') if chrom else 'unlabelled_target_block'
        if out and out[-1]['hi']==start and out[-1]['chrom']==chrom and out[-1]['status']==status:out[-1]['hi']=end
        else:out.append(dict(lo=start,hi=end,chrom=chrom,status=status))
    return out


def labels(job,out,threads):
    catalog=json.loads(Path(job['catalog']).read_text());evidence=json.loads(Path(job['evidence']).read_text())
    regions=defaultdict(list)
    for e in job['events']:regions[e['scaffold']].append((e['lo'],e['hi']))
    observations=defaultdict(list);track_rows=[];hit_rows=[];tiles=[]
    for pair in job['alignments']:
        peer_id=pair['a'] if pair['b']==job['assembly'] else pair['b']
        peer=catalog[peer_id];label=independent_labels(peer,job['sample'])
        print('Projecting local chromosome evidence: '+peer_id,flush=True)
        grouped=local_hits(pair,job,regions)
        for e in job['events']:
            hits=grouped[e['scaffold']]
            projected={mode:local_tracks(hits,label,e['lo'],e['hi'],mode)
                       for mode in ('chromosome','locus_unique')}
            for mode,values in projected.items():
                for b in values:
                    lo=max(e['lo'],b['lo']);hi=min(e['hi'],b['hi'])
                    if hi<=lo:continue
                    track_rows.append(dict(event=e['event'],peer=peer_id,sample=peer['sample'],eligible=peer['eligible'],
                                           mode=mode,lo=lo,hi=hi,chrom=b['chrom'] or '',status=b['status']))
                    if mode=='locus_unique' and b['chrom'] and peer['eligible']:
                        observations[e['event']].append((peer['sample'],lo,hi,b['chrom']))
            for mode,values in projected.items():
                bins=defaultdict(Counter)
                for b in values:
                    if not b['chrom']:continue
                    begin=max(e['lo'],b['lo']);end=min(e['hi'],b['hi'])
                    if end<=begin:continue
                    for i in range((begin-e['lo'])//5000,(end-1-e['lo'])//5000+1):
                        lo=e['lo']+i*5000;hi=min(lo+5000,e['hi'])
                        bins[i][b['chrom']]+=min(hi,end)-max(lo,begin)
                for i,counts in sorted(bins.items()):
                    lo=e['lo']+i*5000;hi=min(lo+5000,e['hi'])
                    for chrom,n in sorted(counts.items()):
                        tiles.append(dict(event=e['event'],peer=peer_id,sample=peer['sample'],mode=mode,
                                          lo=lo,hi=hi,chrom=chrom,aligned_bp=n))
            for f in hits:
                if int(f[2])>=e['hi'] or int(f[3])<=e['lo']:continue
                hit_rows.append(dict(event=e['event'],peer=peer_id,query_lo=f[2],query_hi=f[3],
                                     target=f[5],target_lo=f[7],target_hi=f[8],strand=f[4],
                                     mapq=f[11],identity=int(f[9])/max(1,int(f[10]))))
    islands=[];island_index=defaultdict(list);event_peers=defaultdict(set)
    for r in track_rows:
        event_peers[r['event']].add(r['peer'])
        if r['chrom']:island_index[r['event'],r['peer'],r['mode'],r['chrom']].append((r['lo'],r['hi']))
    for e in job['events']:
        for s in e['sites']:
            for side,chrom,lo,hi in (('right',e['left'],s['point'],e['hi']),
                                    ('left',e['right'],e['lo'],s['point'])):
                for peer in sorted(event_peers[e['event']]):
                    counts={};largest=0;sample=catalog[peer]['sample'];eligible=catalog[peer]['eligible']
                    for mode in ('chromosome','locus_unique'):
                        intervals=union((max(lo,a),min(hi,b)) for a,b in island_index[e['event'],peer,mode,chrom])
                        counts[mode]=sum(b-a for a,b in intervals)
                        if mode=='locus_unique':largest=max((b-a for a,b in intervals),default=0)
                    islands.append(dict(event=e['event'],site=s['id'],peer=peer,sample=sample,eligible=eligible,
                                        side=side,chrom=chrom,chromosome_assigned_bp=counts['chromosome'],
                                        locus_unique_bp=counts['locus_unique'],largest_unique_block_bp=largest))
    masks={e['event']:consensus_masks(observations[e['event']],job['sample']) for e in job['events']}
    save(out/'masks.json',masks)
    write(out/'tracks.tsv',track_rows,['event','peer','sample','eligible','mode','lo','hi','chrom','status'])
    write(out/'chromosome_tiles.tsv',tiles,['event','peer','sample','mode','lo','hi','chrom','aligned_bp'])
    write(out/'alignment_hits.tsv',hit_rows,['event','peer','query_lo','query_hi','target','target_lo','target_hi','strand','mapq','identity'])
    write(out/'opposite_side_islands.tsv',islands,['event','site','peer','sample','eligible','side','chrom','chromosome_assigned_bp','locus_unique_bp','largest_unique_block_bp'])
    save(out/'summary.json',dict(branch='labels',assembly=job['assembly'],
         events={e['event']:dict(unique_consensus_bp=sum(b['hi']-b['lo'] for b in masks[e['event']]),
                                original_graph={k:e0.get(k) for k in ('graph_status','native_left','native_right','native_continuity','direct_native_link')})
                 for e in job['events'] for e0 in evidence['events'] if e0['id']==e['event']},
         limitations=['Unique means unique among emitted assembly alignments, not proven genome-wide single copy.',
                      'No repeat library annotation is inferred from absent alignments.',
                      'Graph context is retained from the original native-path assay; it is not independent read evidence.']))


def sam_records(path):
    with open(path) as src:
        for line in src:
            if not line.startswith('@'):yield line.rstrip().split('\t')


def span_record(f,point,anchor):
    return any(lo<=point-anchor and hi>=point+anchor for lo,hi in molecule_chains(f))


def tag(f,name,default=0):
    return next((int(t.split(':',2)[2]) for t in f[11:] if t.startswith(name+':i:')),default)


def alignment_summary(f):
    ops=re.findall(r'(\d+)([MIDNSHP=X])',f[5])
    aligned=sum(int(n) for n,op in ops if op in 'MI=X')
    return dict(score=tag(f,'AS'),identity=1-tag(f,'NM')/max(1,aligned),
                aligned=aligned,clipped=sum(int(n) for n,op in ops if op in 'SH'),
                large_indels=sum(int(n)>50 for n,op in ops if op in 'IDN'))


def arm_coverage(f,point,mask,left,right):
    """Unique chromosome-specific bases on each side within one uninterrupted chain."""
    results=[]
    for lo,hi in molecule_chains(f):
        if lo>point or hi<point:continue
        l=overlap([(b['lo'],b['hi']) for b in mask if b['chrom']==left],lo,point)
        r=overlap([(b['lo'],b['hi']) for b in mask if b['chrom']==right],point,hi)
        results.append((l,r))
    return max(results,key=lambda x:min(x),default=(0,0))


def minimum_anchor_span(mask,point,left,right,anchor=1000):
    """Distance to accumulate anchor unique bases on each side, allowing islands."""
    ls=union((b['lo'],min(point,b['hi'])) for b in mask if b['chrom']==left and b['lo']<point)
    rs=union((max(point,b['lo']),b['hi']) for b in mask if b['chrom']==right and b['hi']>point)
    need=anchor;start=None;end=None
    for lo,hi in reversed(ls):
        if hi-lo>=need:start=hi-need;break
        need-=hi-lo
    need=anchor
    for lo,hi in rs:
        if hi-lo>=need:end=lo+need;break
        need-=hi-lo
    return end-start if start is not None and end is not None else None


def hifi(job,out,threads,masks_path):
    # Merge overlapping extraction windows: each molecule enters competitive mapping once.
    masks=json.loads(Path(masks_path).read_text());regions=defaultdict(list)
    for e in job['events']:regions[e['scaffold']].append((max(0,e['lo']-100000),min(job['sizes'][e['scaffold']],e['hi']+100000)))
    bed=out/'regions.bed'
    bed.write_text(''.join(f'{s}\t{lo}\t{hi}\n' for s,rs in regions.items() for lo,hi in union(rs)))
    sam=out/'original.sam';run(['samtools','view','-h','-M','-L',bed,job['bam']],sam)
    reads={};records=[];profiles=[];site_rows=[]
    for f in sam_records(sam):
        if int(f[1])&(4|256|512|1024|2048):continue
        records.append(f)
        if f[9]!='*':
            seq=f[9]
            if int(f[1])&16:seq=seq.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1]
            if f[0] in reads and reads[f[0]]!=seq:raise ValueError('Inconsistent primary molecule sequences')
            reads[f[0]]=seq
    fasta=out/'reads.fa'
    with fasta.open('w') as dst:
        for name,seq in reads.items():dst.write(f'>{name}\n{seq}\n')
    for e in job['events']:
        points=list(range(e['lo'],e['hi'],1000));counts={k:[0]*(len(points)+1) for k in (0,20,30)}
        for f in records:
            if f[2]!=e['scaffold']:continue
            copy=list(f);copy[4]=str(max(30,int(f[4])))
            if not qualified_molecule(copy):continue
            # Merge chain-supported probe intervals before counting this molecule.
            ranges=union((bisect_left(points,lo+1000),bisect_right(points,hi-1000)) for lo,hi in molecule_chains(f))
            for tier in counts:
                if int(f[4])<tier:continue
                for lo,hi in ranges:counts[tier][lo]+=1;counts[tier][hi]-=1
        totals={tier:0 for tier in counts}
        for i,p in enumerate(points):
            for tier in totals:totals[tier]+=counts[tier][i]
            profiles.append(dict(event=e['event'],position=p,**{'mapq'+str(t):totals[t] for t in totals}))
        for s in e['sites']:
            for anchor in (1000,5000):
                bridgers=[]
                for f in records:
                    copy=list(f);copy[4]=str(max(30,int(f[4])))
                    if f[2]==e['scaffold'] and qualified_molecule(copy) and span_record(f,s['point'],anchor):
                        bridgers.append(f)
                site_rows.append(dict(event=e['event'],site=s['id'],point=s['point'],anchor_bp=anchor,
                                      mapq0=len({f[0] for f in bridgers}),
                                      mapq20=len({f[0] for f in bridgers if int(f[4])>=20}),
                                      mapq30=len({f[0] for f in bridgers if int(f[4])>=30})))
    write(out/'spanning_profile.tsv',profiles,['event','position','mapq0','mapq20','mapq30'])
    write(out/'original_bridges.tsv',site_rows,['event','site','point','anchor_bp','mapq0','mapq20','mapq30'])
    # Separate indexes preserve each assembly's placement competition. A homologous
    # placement in the other haplotype is reported, not used as an automatic veto.
    candidate_reads=defaultdict(list)
    for target_assembly,target in job['competitive_fastas'].items():
        mapped=out/(target_assembly+'.sam')
        if reads:
            run(['minimap2','-ax','map-hifi','-I','8G','--secondary=yes','-N','50','-p','0.5',
                 '-t',threads,target,fasta],mapped)
        else:mapped.write_text('')
        for f in sam_records(mapped):
            if int(f[1])&4:continue
            info=alignment_summary(f)
            if info['identity']<.9:continue
            candidate_reads[f[0]].append((target_assembly,f,info))
    bridges=[];placements=[];reach=[]
    for e in job['events']:
        for s in e['sites']:
            mask=masks[e['event']]
            required_span=minimum_anchor_span(mask,s['point'],e['left'],e['right'])
            local_lengths=sorted(len(reads[f[0]]) for f in records if f[0] in reads and f[2]==e['scaffold'] and any(lo<s['point']+100000 and hi>s['point']-100000 for lo,hi in molecule_chains(f)))
            reach.append(dict(event=e['event'],site=s['id'],minimum_unique_anchor_span_bp=required_span,
                              local_read_count=len(local_lengths),local_read_median=statistics.median(local_lengths) if local_lengths else None,
                              local_read_p95=local_lengths[min(len(local_lengths)-1,int(.95*len(local_lengths)))] if local_lengths else None,
                              local_reads_long_enough=sum(n>=required_span for n in local_lengths) if required_span is not None else None))
            for name,rs in candidate_reads.items():
                for target,f,info in rs:
                    if target!=job['assembly'] or f[2]!=e['scaffold'] or not span_record(f,s['point'],1000):continue
                    if int(f[1])&(256|2048):continue
                    alternatives=[(ta,g,gi) for ta,g,gi in rs if ta==target and
                                  (g[2]!=f[2] or bool(int(g[1])&16)!=bool(int(f[1])&16) or abs(int(g[3])-int(f[3]))>1000)
                                  and gi['aligned']>=.9*info['aligned'] and gi['score']>=.95*info['score']]
                    homologs=[(ta,g,gi) for ta,g,gi in rs if ta!=target and gi['aligned']>=.9*info['aligned'] and gi['score']>=.95*info['score']]
                    left_bp,right_bp=arm_coverage(f,s['point'],mask,e['left'],e['right'])
                    bridges.append(dict(event=e['event'],site=s['id'],point=s['point'],read=name,mapq=int(f[4]),
                                        **info,anchor5kb=span_record(f,s['point'],5000),
                                        unique_left_bp=left_bp,unique_right_bp=right_bp,
                                        same_haplotype_alternatives=len(alternatives),other_haplotype_alternatives=len(homologs)))
                    for ta,g,gi in rs:
                        placements.append(dict(event=e['event'],site=s['id'],read=name,assembly=ta,
                                               scaffold=g[2],start=int(g[3])-1,flag=int(g[1]),mapq=int(g[4]),cigar=g[5],**gi))
    write(out/'competitive_bridges.tsv',bridges,['event','site','point','read','mapq','score','identity','aligned','clipped','large_indels','anchor5kb','unique_left_bp','unique_right_bp','same_haplotype_alternatives','other_haplotype_alternatives'])
    write(out/'unique_anchor_reach.tsv',reach,['event','site','minimum_unique_anchor_span_bp','local_read_count','local_read_median','local_read_p95','local_reads_long_enough'])
    write(out/'bridge_placements.tsv',placements,['event','site','read','assembly','scaffold','start','flag','mapq','cigar','score','identity','aligned','clipped','large_indels'])
    save(out/'summary.json',dict(branch='hifi',assembly=job['assembly'],selected_molecules=len(reads),
         branch_status='measured',limitations=['Competing placements are bounded by minimap2 emitted alternatives.',
         'Other-haplotype placement alone does not invalidate a bridge.',
         'Exact repeat-spanning reach is assessed against unique-flank distances in the report.',
         'Read selection uses primary placements overlapping these windows; reads with only secondary local placements are not exhaustively recovered.']))
    # Retain these large, diagnostic-only products in the work directory, not the packet.
    for p in out.glob('*.sam'):p.unlink()
    fasta.unlink();bed.unlink()


def library_id(qname,libraries):
    found={r['library_id'] for r in libraries if qname.startswith(r['qname_prefix'])}
    if len(found)!=1:raise ValueError('Missing or ambiguous Hi-C library prefix')
    return next(iter(found))


def window_counts(matrix,point,flank):
    """Symmetric counts with diagonal included once; point and flank are bin indices."""
    import numpy as np
    if point-flank<0 or point+flank>len(matrix):return None
    left=matrix[point-flank:point,point-flank:point]
    right=matrix[point:point+flank,point:point+flank]
    l=int((left.sum(dtype=np.uint64)+np.trace(left,dtype=np.uint64))//2)
    r=int((right.sum(dtype=np.uint64)+np.trace(right,dtype=np.uint64))//2)
    cross=int(matrix[point-flank:point,point:point+flank].sum(dtype=np.uint64))
    return dict(cross=cross,left_within=l,right_within=r,
                ratio=cross/math.sqrt(l*r) if l and r else None)


def aggregate_matrix(base,factor):
    """Coarsen a symmetric pair matrix without counting within-bin pairs twice."""
    import numpy as np
    n=len(base);extra=(-n)%factor
    padded=np.pad(base,((0,extra),(0,extra))) if extra else base
    coarse=padded.reshape(len(padded)//factor,factor,len(padded)//factor,factor).sum(axis=(1,3),dtype=np.uint64)
    for i in range(len(coarse)):
        block=padded[i*factor:(i+1)*factor,i*factor:(i+1)*factor]
        coarse[i,i]=(block.sum(dtype=np.uint64)+np.trace(block,dtype=np.uint64))//2
    return coarse


def mask_bins(mask,lo,n,resolution,chrom,minimum=.8):
    return [i for i in range(n) if overlap([(b['lo'],b['hi']) for b in mask if b['chrom']==chrom],
                                         lo+i*resolution,lo+(i+1)*resolution)>=minimum*resolution]


def hic(job,out,threads,masks_path):
    import numpy as np
    masks=json.loads(Path(masks_path).read_text());libraries=table(job['paths']['libraries'])
    prefixes=[r['qname_prefix'] for r in libraries]
    if not prefixes or any(not p for p in prefixes) or any(p.startswith(q) for i,p in enumerate(prefixes) for j,q in enumerate(prefixes) if i!=j):
        raise ValueError('Missing or overlapping library prefixes')
    libs=sorted({r['library_id'] for r in libraries})
    _,placements,lengths=agp(job['paths']['agp'])
    for sc,size in lengths.items():
        if sc not in job['sizes'] or size!=job['sizes'][sc]:raise ValueError('AGP/assessment lengths differ')
    def locate(sc,pos):
        found=[v for v in placements.get(sc,[]) if v[3]<=pos<v[4]]
        if len(found)!=1:return None
        v=found[0];return v[0],v[1]+pos-v[3] if v[5]=='+' else v[2]-1-(pos-v[3])
    # One smallest-resolution matrix per event and library; aggregate exactly later.
    resolution=5000;matrices={};geometry={};by_scaffold=defaultdict(list)
    gaps=literal_gaps(job['paths']['fasta'],minimum=1)
    for e in job['events']:
        lo=e['lo']//100000*100000;hi=min(job['sizes'][e['scaffold']],math.ceil(e['hi']/100000)*100000)
        n=math.ceil((hi-lo)/resolution);geometry[e['event']]=(lo,hi,n)
        by_scaffold[e['scaffold']].append(e)
        for lib in libs:matrices[e['event'],lib]=np.zeros((n,n),dtype=np.uint32)
    audit=Counter();opener=gzip.open if job['paths']['pairs'].endswith('.gz') else open
    with opener(job['paths']['pairs'],'rt') as src:
        for line in src:
            if not line.strip() or line.startswith('#'):continue
            f=line.rstrip().split('\t');audit['input_pairs']+=1
            if len(f)<8:raise ValueError('Malformed pairs row')
            if f[7]!='UU':continue
            lib=library_id(f[0],libraries);x=locate(f[1],int(f[2])-1);y=locate(f[3],int(f[4])-1)
            if x is None or y is None:audit['unliftable']+=1;continue
            audit['liftable']+=1
            if x[0]!=y[0]:continue
            for e in by_scaffold.get(x[0],[]):
                lo,hi,n=geometry[e['event']]
                if lo<=x[1]<hi and lo<=y[1]<hi:
                    i=(x[1]-lo)//resolution;j=(y[1]-lo)//resolution
                    m=matrices[e['event'],lib];m[i,j]+=1
                    if i!=j:m[j,i]+=1
    if audit['input_pairs']!=job['hic_expected']:raise ValueError('Pairs row count differs from source evidence')
    profiles=[];anchors=[];controls=[];marginals=[];map_records=[]
    for e in job['events']:
        lo,hi,n=geometry[e['event']]
        mask=masks[e['event']]
        unique_bins={chrom:set(mask_bins(mask,lo,n,resolution,chrom)) for chrom in (e['left'],e['right'])}
        for lib in libs+['POOLED']:
            base=sum((matrices[e['event'],l].astype(np.uint64) for l in libs),np.zeros((n,n),dtype=np.uint64)) if lib=='POOLED' else matrices[e['event'],lib]
            for res in (10000,25000,50000,100000):
                factor=res//resolution;coarse=aggregate_matrix(base,factor)
                filename=f"{e['event']}.{lib}.{res}.npz";np.savez_compressed(out/filename,matrix=coarse,lo=lo,resolution=res)
                map_records.append(dict(event=e['event'],library=lib,resolution=res,file=filename,lo=lo,hi=hi))
            for i in range(n):
                marginals.append(dict(event=e['event'],library=lib,lo=lo+i*resolution,hi=min(hi,lo+(i+1)*resolution),
                                      local_contact_ends=int(base[i].sum(dtype=np.uint64)+base[i,i])))
            for flank_bp in (100000,300000,500000):
                flank=flank_bp//resolution
                for i in range(flank,n-flank+1):
                    if i%2:continue  # 10 kb scan steps, fixed flank sizes.
                    point=lo+i*resolution
                    if point+flank_bp>hi:continue
                    c=window_counts(base,i,flank)
                    profiles.append(dict(event=e['event'],library=lib,point=point,flank_bp=flank_bp,**c))
                    # Controls are gap-free in the retained reference and have unique,
                    # same-chromosome cohort support across both counting windows.
                    same=[chrom for chrom in (e['left'],e['right']) if
                          overlap([(b['lo'],b['hi']) for b in mask if b['chrom']==chrom],
                                  point-flank_bp,point+flank_bp)>=1.6*flank_bp]
                    gap_free=not any(g['scaffold']==e['scaffold'] and g['lo']<point+flank_bp and g['hi']>point-flank_bp for g in gaps)
                    if same and gap_free and all(abs(point-s['point'])>2*flank_bp for s in e['sites']):
                        controls.append(dict(event=e['event'],library=lib,point=point,flank_bp=flank_bp,chrom=same[0],**c))
            # Unique arm windows: fixed sizes, nearest usable 10 kb bins on the proper
            # side of each site. Their actual coordinates and separation are retained.
            for s in e['sites']:
                left=[i for i in sorted(unique_bins[e['left']]) if lo+(i+1)*resolution<=s['point']]
                right=[i for i in sorted(unique_bins[e['right']]) if lo+i*resolution>=s['point']]
                for width in (50000,100000,300000):
                    need=width//resolution;ls=left[-need:];rs=right[:need]
                    row=dict(event=e['event'],site=s['id'],point=s['point'],library=lib,anchor_bp=width,
                             left_bins=len(ls),right_bins=len(rs),status='insufficient_unique_arm_sequence')
                    if len(ls)==len(rs)==need:
                        cross=int(base[np.ix_(ls,rs)].sum(dtype=np.uint64))
                        lmat=base[np.ix_(ls,ls)];rmat=base[np.ix_(rs,rs)]
                        lw=int((lmat.sum(dtype=np.uint64)+np.trace(lmat,dtype=np.uint64))//2)
                        rw=int((rmat.sum(dtype=np.uint64)+np.trace(rmat,dtype=np.uint64))//2)
                        row.update(status='measured',left_lo=lo+min(ls)*resolution,left_hi=lo+(max(ls)+1)*resolution,
                                   right_lo=lo+min(rs)*resolution,right_hi=lo+(max(rs)+1)*resolution,
                                   cross=cross,left_within=lw,right_within=rw,ratio=cross/math.sqrt(lw*rw) if lw and rw else None)
                        row['left_intervals']=';'.join(f'{a}:{b}' for a,b in union((lo+i*resolution,lo+(i+1)*resolution) for i in ls))
                        row['right_intervals']=';'.join(f'{a}:{b}' for a,b in union((lo+i*resolution,lo+(i+1)*resolution) for i in rs))
                        comparisons=[]
                        for shift in range(-n,n,max(1,width//resolution)):
                            if shift==0:continue
                            lshift=[i+shift for i in ls];rshift=[i+shift for i in rs]
                            if min(lshift)<0 or max(rshift)>=n:continue
                            center=lo+(s['point']-lo)+shift*resolution
                            if any(abs(center-t['point'])<2*width for t in e['sites']):continue
                            selected=set(lshift+rshift)
                            chroms=[c for c in (e['left'],e['right']) if selected.issubset(unique_bins[c])]
                            if not chroms:continue
                            lm=base[np.ix_(lshift,lshift)];rm=base[np.ix_(rshift,rshift)]
                            cl=int((lm.sum(dtype=np.uint64)+np.trace(lm,dtype=np.uint64))//2)
                            cr=int((rm.sum(dtype=np.uint64)+np.trace(rm,dtype=np.uint64))//2)
                            cx=int(base[np.ix_(lshift,rshift)].sum(dtype=np.uint64))
                            if min(cl,cr,lw,rw)<100 or not (.1<=cl/lw<=10 and .1<=cr/rw<=10):continue
                            value=cx/math.sqrt(cl*cr);comparisons.append(value)
                            controls.append(dict(event=e['event'],site=s['id'],library=lib,point=center,flank_bp=width,
                                                 chrom=chroms[0],cross=cx,left_within=cl,right_within=cr,ratio=value,kind='translated_unique_arms'))
                        expected=statistics.median(comparisons)*math.sqrt(lw*rw) if len(comparisons)>=5 else None
                        row.update(matched_controls=len(comparisons),expected_cross=expected,
                                   observed_expected=cross/expected if expected else None)
                    anchors.append(row)
    # Replay the old 100 kb / 300 kb-flank assay exactly wherever the local
    # window includes the full geometry. This catches lift and binning mistakes.
    lookup={(r['event'],r['library'],r['point']):r for r in profiles if r['flank_bp']==300000}
    old=json.loads(Path(job['evidence']).read_text());replayed=0
    for e in old['events']:
        for lib,rows in e.get('contact_profiles',{}).items():
            for row in rows:
                current=lookup.get((e['id'],lib,row['position']))
                if current is None:continue
                for field in ('cross','left_within','right_within'):
                    if current[field]!=row[field]:raise ValueError(f'100 kb replay differs: {e["id"]} {lib} {row["position"]} {field}')
                replayed+=1
    if not replayed:raise ValueError('No original contact measurements could be replayed')
    write(out/'contact_profile.tsv',profiles,['event','library','point','flank_bp','cross','left_within','right_within','ratio'])
    write(out/'unique_arm_contacts.tsv',anchors,['event','site','point','library','anchor_bp','status','left_bins','right_bins','left_lo','left_hi','right_lo','right_hi','left_intervals','right_intervals','cross','left_within','right_within','ratio','matched_controls','expected_cross','observed_expected'])
    write(out/'controls.tsv',controls,['event','site','kind','library','point','flank_bp','chrom','cross','left_within','right_within','ratio'])
    write(out/'marginals.tsv',marginals,['event','library','lo','hi','local_contact_ends'])
    save(out/'summary.json',dict(branch='hic',assembly=job['assembly'],audit=dict(audit),maps=map_records,replayed_measurements=replayed,
         limitations=['Local marginals count contacts retained within the displayed window, not genome-wide coverage.',
         'Pooled measurements supplement separate libraries and are not an independent replicate.',
         'Overlapping scan positions and control windows are correlated; no independent-bin significance is claimed.']))


def shade_masks(ax,mask,colors,alpha):
    """One collection per chromosome, never one plot object per alignment block."""
    from matplotlib.collections import PolyCollection
    grouped=defaultdict(list)
    for b in mask:
        lo=b['lo']/1e6;hi=b['hi']/1e6
        grouped[colors.get(b['chrom'],'grey')].append([(lo,0),(hi,0),(hi,1),(lo,1)])
    for color,vertices in grouped.items():
        ax.add_collection(PolyCollection(vertices,facecolors=color,edgecolors='none',alpha=alpha,
                                        transform=ax.get_xaxis_transform()),autolim=False)


def report(a):
    """One review packet, with descriptive evidence and explicit unselected actions."""
    import shutil
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    manifest=json.loads(a.manifest.read_text());a.out.mkdir(exist_ok=False)
    lines=['# Four-join investigation','',
           'Read-only evidence review. **No cuts are selected or applied.**', '',
           'Coordinates below refer to the checksum-verified pre-finishing assessment FASTA. '
           'Final chromosome names are context only: finishing can shift coordinates, so they are not an exact coordinate lift.', '',
           'Blue and orange identify the two chromosome arms within each case. '
           'HiFi spanning molecules are long reads; Hi-C contacts are pairs of short reads.', '',
           '## Cases','', '| Assembly / case | Chromosomes | Role | Evidence |', '|---|---|---|---|']
    case_lines=[];actions=[];summaries=[]
    for job in manifest['jobs']:
        assembly=job['assembly'];dest=a.out/assembly;dest.mkdir()
        for branch in ('labels','hifi','hic'):
            folder=next(p for p in a.packets if p.name==assembly+'.'+branch)
            shutil.copytree(folder,dest/branch)
        hifi=table(dest/'hifi'/'spanning_profile.tsv');original=table(dest/'hifi'/'original_bridges.tsv')
        bridges=table(dest/'hifi'/'competitive_bridges.tsv');contacts=table(dest/'hic'/'contact_profile.tsv')
        reach=table(dest/'hifi'/'unique_anchor_reach.tsv')
        tracks_data=table(dest/'labels'/'tracks.tsv');anchors=table(dest/'hic'/'unique_arm_contacts.tsv')
        masks=json.loads((dest/'labels'/'masks.json').read_text());maps=json.loads((dest/'hic'/'summary.json').read_text())['maps']
        names=table(job['final_map']);write(dest/'current_names.tsv',names,list(names[0]))
        for e in job['events']:
            key=assembly+'-'+e['event'];anchor=key.lower().replace('_','-')
            role='Target' if e['role']=='target' else 'Previously reviewed break (comparison)'
            lines.append(f"| {assembly} / {e['event']} | {e['left']} → {e['right']} | {role} | [Review](#{anchor}) |")
            case_lines += ['',f'<a id="{anchor}"></a>',f"## {assembly} / {e['event']}: {e['left']} → {e['right']}",'',
                           '### What this batch can establish','',
                           'Compare the chromosome transition with every measured support trough. '
                           'A small alignment island beyond a trough does not veto a cut; inspect its unique-locus support. '
                           'A trough away from the arm transition may instead be a separate local assembly problem.', '',
                           f"Assessment scaffold: `{e['scaffold']}`; transition interval **{e['transition_lo']:,}–{e['transition_hi']:,} bp**.", '']
            matches=[r['new_name'] for r in names if r['old_name']==e['scaffold'] or r['old_name'].startswith(e['scaffold']+'_sub_')]
            case_lines += ['Current named descendants: '+(', '.join('`'+v+'`' for v in matches) or 'No matching descendant in supplied name map.')+'.','',
                           '### Sites and molecule evidence','',
                           '| Site | Position, bp | Original HiFi ≥20 / ≥30 | HiFi spans after checking alternatives | Spans joining ≥1 kb unique sequence from both chromosomes |',
                           '|---|---:|---:|---:|---:|']
            for s in e['sites']:
                old=next(r for r in original if r['site']==s['id'] and r['anchor_bp']=='1000')
                credible=[r for r in bridges if r['site']==s['id'] and int(r['mapq'])>=30 and int(r['same_haplotype_alternatives'])==0]
                n=len({r['read'] for r in credible});unique=len({r['read'] for r in credible if min(int(r['unique_left_bp']),int(r['unique_right_bp']))>=1000})
                kind={'transition_midpoint':'transition midpoint','regional_trough':'support trough','gap_option':'scaffold gap'}[s['kind']]
                case_lines.append(f"| [{s['id']}]({assembly}/{e['event']}.{s['id']}.detail.png) ({kind}) | {s['point']:,} | {old['mapq20']} / {old['mapq30']} | {n} | {unique} |")
                actions.append(dict(assembly=assembly,id=s['id'],scaffold=e['scaffold'],selected='NO',action='BREAK_PROBABLE_MISJOIN',
                                    cut_bp='',assay_probe_bp=s['point'],reviewer='',reason='',coordinate_stage='pre_finishing',
                                    assessment_sha256=job['assessment_sha256'],decision_source='review',evidence_packet_id=key,
                                    localization_status='localized',auto_eligible='no'))
                summaries.append(dict(assembly=assembly,event=e['event'],site=s['id'],role=e['role'],position=s['point'],
                                      original_mapq20=old['mapq20'],original_mapq30=old['mapq30'],competitive_credible=n,unique_arm_bridges=unique))
            case_lines += ['', 'The transition midpoint is a measurement location, **not a proposed cut**. '
                           'Competitive counts are descriptive, not a haplotype phasing test. '
                           'Other-haplotype similarity alone does not invalidate a bridge. Click a site ID for its detailed view.', '',
                           '### Can a read reach the unique chromosome anchors?','',
                           '| Site | Minimum span needed | Local read median / 95th percentile | Local reads long enough |',
                           '|---|---:|---:|---:|']
            for r in reach:
                if r['event']!=e['event']:continue
                case_lines.append(f"| {r['site']} | {r['minimum_unique_anchor_span_bp'] or 'No usable anchors'} | {r['local_read_median']} / {r['local_read_p95']} | {r['local_reads_long_enough'] or '0'} |")
            case_lines += ['', 'Read length feasibility is not an expected bridge count: start positions, coverage and phasing also matter.', '',
                           '### Chromosome-labelled islands beyond each site','',
                           '| Site | Chromosome on opposite side | Independent individuals with unique support | Largest unique block in a peer |',
                           '|---|---|---:|---:|']
            island_rows=table(dest/'labels'/'opposite_side_islands.tsv')
            for s in e['sites']:
                for chrom,side in ((e['left'],'right'),(e['right'],'left')):
                    rs=[r for r in island_rows if r['event']==e['event'] and r['site']==s['id'] and r['chrom']==chrom and r['side']==side]
                    independent={r['sample'] for r in rs if r['sample'] and r['sample']!=job['sample'] and r['eligible']=='True' and int(r['locus_unique_bp'])>0}
                    case_lines.append(f"| {s['id']} | {chrom}, {side} side | {len(independent)} | {max((int(r['largest_unique_block_bp']) for r in rs),default=0):,} bp |")
            case_lines += ['', 'These are observations, not a rule that a cut must follow the last alignment. '
                           'A long unique block needs a different explanation from a short multiply placed island. '
                           f'[All island lengths and per-peer evidence]({assembly}/labels/opposite_side_islands.tsv).','',
                           '### Aligned evidence','',f'![Evidence]({assembly}/{e["event"]}.evidence.png)','']
            fig,axes=plt.subplots(4,1,figsize=(15,12),sharex=True,gridspec_kw={'height_ratios':[2,1,1,1]},layout='constrained')
            peers=sorted({r['peer'] for r in tracks_data if r['event']==e['event']})
            ordered_chroms=sorted((e['left'],e['right']),key=lambda c:[int(v) if v.isdigit() else v for v in re.split(r'(\d+)',c)])
            colors=dict(zip(ordered_chroms,('#0072B2','#D55E00')))
            for i,peer in enumerate(peers):
                grouped=defaultdict(list)
                for r in tracks_data:
                    if r['event']!=e['event'] or r['peer']!=peer or r['mode']!='locus_unique' or not r['chrom']:continue
                    grouped[colors.get(r['chrom'],'#999999')].append((int(r['lo'])/1e6,(int(r['hi'])-int(r['lo']))/1e6))
                for color,intervals in grouped.items():axes[0].broken_barh(intervals,(i-.35,.7),facecolors=color)
            axes[0].set_yticks(range(len(peers)),peers);axes[0].set_title('Unique-locus chromosome assignments; blank = unresolved, not absent sequence',pad=45,fontsize=11)
            rs=[r for r in hifi if r['event']==e['event']]
            for tier,c in (('mapq20','#0072B2'),('mapq30','#D55E00')):
                axes[1].plot([int(r['position'])/1e6 for r in rs],[int(r[tier]) for r in rs],label=tier,color=c,lw=1)
            axes[1].set_ylabel('HiFi spanning\nmolecules');axes[1].legend()
            for lib in sorted({r['library'] for r in contacts}-{'POOLED'}):
                rs=[r for r in contacts if r['event']==e['event'] and r['library']==lib and r['flank_bp']=='300000']
                x=[int(r['point'])/1e6 for r in rs]
                y=[float(r['ratio']) if r['ratio'] else math.nan for r in rs]
                axes[2].plot(x,y,label=lib,lw=1)
                low=[i for i,r in enumerate(rs) if min(int(r['left_within']),int(r['right_within']))<100]
                axes[2].scatter([x[i] for i in low],[y[i] for i in low],marker='x',color='grey',s=12)
                axes[3].plot(x,[min(int(r['left_within']),int(r['right_within'])) for r in rs],label=lib,lw=1)
            axes[2].set_ylabel('Hi-C cross /\n√(within L × R)');axes[2].legend();axes[2].set_ylim(bottom=0)
            axes[3].set_ylabel('Weaker flank\nwithin-pairs');axes[3].set_yscale('symlog',linthresh=1);axes[3].set_ylim(bottom=0)
            for ax in axes:
                if ax is not axes[0]:
                    shade_masks(ax,masks[e['event']],colors,.08)
                ax.axvspan(e['transition_lo']/1e6,e['transition_hi']/1e6,color='#E69F00',alpha=.12)
                for s in e['sites']:ax.axvline(s['point']/1e6,color='#882255',ls='--',lw=.7,alpha=.7)
                ax.set_facecolor('white');ax.grid(alpha=.15)
            for i,s in enumerate(e['sites']):axes[0].text(s['point']/1e6,1.02+(i%2)*.07,s['id'],transform=axes[0].get_xaxis_transform(),ha='center',fontsize=7)
            axes[-1].set_xlabel('Assessment scaffold coordinate (Mb)');axes[-1].set_xlim(e['lo']/1e6,e['hi']/1e6)
            fig.suptitle(f"{assembly} / {e['event']}: {e['left']} → {e['right']} · blue={ordered_chroms[0]}, orange={ordered_chroms[1]}")
            fig.savefig(dest/(e['event']+'.evidence.png'),dpi=140,facecolor='white');plt.close(fig)
            # Every measured site has a detailed view, even if it is outside the
            # chromosome-transition interval or has no recommended exact cut.
            for s in e['sites']:
                fig,axs=plt.subplots(3,1,figsize=(12,7),sharex=True,layout='constrained')
                for tier,c in (('mapq20','#0072B2'),('mapq30','#D55E00')):
                    rs=[r for r in hifi if r['event']==e['event']]
                    axs[0].plot([int(r['position'])/1e6 for r in rs],[int(r[tier]) for r in rs],color=c,label=tier,lw=1)
                axs[0].legend();axs[0].set_ylabel('HiFi molecules')
                for lib in sorted({r['library'] for r in contacts}-{'POOLED'}):
                    rs=[r for r in contacts if r['event']==e['event'] and r['library']==lib and r['flank_bp']=='100000']
                    axs[1].plot([int(r['point'])/1e6 for r in rs],[float(r['ratio']) if r['ratio'] else math.nan for r in rs],label=lib,lw=1)
                    axs[2].plot([int(r['point'])/1e6 for r in rs],[min(int(r['left_within']),int(r['right_within'])) for r in rs],label=lib,lw=1)
                axs[1].set_ylabel('Hi-C ratio\n100 kb flanks');axs[1].legend();axs[1].set_ylim(bottom=0);axs[2].set_ylabel('Weaker flank\nwithin-pairs')
                for ax in axs:
                    shade_masks(ax,masks[e['event']],colors,.12)
                    ax.axvline(s['point']/1e6,color='#882255',ls='--');ax.grid(alpha=.15);ax.set_facecolor('white')
                axs[-1].set_xlim(max(e['lo'],s['point']-100000)/1e6,min(e['hi'],s['point']+100000)/1e6)
                axs[-1].set_xlabel('Assessment scaffold coordinate (Mb)');fig.suptitle(f"{assembly} / {s['id']} · probe {s['point']:,} bp")
                fig.savefig(dest/f"{e['event']}.{s['id']}.detail.png",dpi=130,facecolor='white');plt.close(fig)
            case_lines += ['Grey crosses mark contact ratios with fewer than 100 within-flank pairs on either side; '
                           'this is an observability flag, not a significance threshold. '
                           'The fourth panel shows the actual weaker-flank count. Dashed lines mark every tested site.', '',
                           '### Unique chromosome-arm contacts','',
                           '| Site | Library | Unique sequence per arm | Cross pairs | Within left / right | Expected cross pairs from matched controls |',
                           '|---|---|---:|---:|---:|---:|']
            for r in anchors:
                if r['event']!=e['event'] or r['anchor_bp']!='100000' or r['library']=='POOLED':continue
                if r['status']!='measured':case_lines.append(f"| {r['site']} | {r['library']} | Insufficient unique sequence | — | — | — |")
                else:case_lines.append(f"| {r['site']} | {r['library']} | 100 kb | {r['cross']} | {r['left_within']} / {r['right_within']} | {round(float(r['expected_cross']),1) if r['expected_cross'] else 'Insufficient matched controls'} |")
            case_lines += ['', 'Selected unique bins can be discontinuous. Full coordinates and other anchor sizes are in '
                           f'[{assembly}/hic/unique_arm_contacts.tsv]({assembly}/hic/unique_arm_contacts.tsv). '
                           'Longer separation naturally reduces contacts; raw cross counts alone do not rank cuts.', '',
                           '### Local contact maps','']
            for record in maps:
                if record['event']!=e['event'] or record['library']=='POOLED':continue
                arr=np.load(dest/'hic'/record['file']);m=arr['matrix'];lo=record['lo'];res=record['resolution'];extent=[lo/1e6,(lo+len(m)*res)/1e6]*2
                fig,ax=plt.subplots(figsize=(7,6),layout='constrained')
                im=ax.imshow(np.log1p(m),origin='lower',extent=extent,cmap='magma',interpolation='nearest')
                for s in e['sites']:
                    ax.axvline(s['point']/1e6,color='cyan',lw=.5);ax.axhline(s['point']/1e6,color='cyan',lw=.5)
                ax.set_title(f"{e['event']} · {record['library']} · {res//1000} kb · raw counts")
                ax.set_xlabel('Assessment coordinate (Mb)');ax.set_ylabel('Assessment coordinate (Mb)');fig.colorbar(im,ax=ax,label='log(1 + pairs)')
                name=f"{e['event']}.{record['library']}.{res}.png";fig.savefig(dest/name,dpi=120,facecolor='white');plt.close(fig)
                if res==25000:case_lines += [f'![Contact map]({assembly}/{name})','']
                else:case_lines += [f'[{record["library"]}: {res//1000} kb contact map]({assembly}/{name})','']
            old_figure=Path(job['evidence']).parent/'plots'/(e['event']+'.png')
            if old_figure.exists():
                shutil.copyfile(old_figure,dest/(e['event']+'.original.png'))
                case_lines += ['<details>','<summary>Original overview, including depth and telomere evidence</summary>','',
                               f'![Original evidence overview]({assembly}/{e["event"]}.original.png)','','</details>','']
            case_lines += ['### Interpretation and limitations','',
                           'A break is favored by a chromosome-arm change plus loss of credible molecule continuity and '
                           'a contact deficit where both arms are observable. Retention is favored by multiple credible '
                           'bridging molecules with adequate anchors and no comparably good same-haplotype alternatives. '
                           'A trough entirely within one uniquely assigned arm may be a separate local problem.', '',
                           '**Decision remains for review. This diagnostic does not apply a new production scoring policy.**', '',
                           f'- [Molecule placements]({assembly}/hifi/bridge_placements.tsv)',
                           f'- [Chromosome evidence by 5 kb tile]({assembly}/labels/chromosome_tiles.tsv)',
                           f'- [Raw alignment locations and orientation]({assembly}/labels/alignment_hits.tsv)',
                           f'- [Contact scan]({assembly}/hic/contact_profile.tsv)',
                           f'- [Same-chromosome control measurements]({assembly}/hic/controls.tsv)',
                           f'- [Native graph context and limitations]({assembly}/labels/summary.json)','']
    lines += case_lines
    lines += ['## Manual cut worksheet','',
              '[manual-cut-worksheet.tsv](manual-cut-worksheet.tsv) contains **unselected** rows with blank cut coordinates. '
              'Assay probes are not preselected cuts. After review, enter the chosen exact coordinate, reviewer and reason. '
              'Use the production manual-review interface to validate and apply selected actions; this investigation never edits FASTA.', '',
              '## Scope','',
              'This packet audits alignment ambiguity using retained assembly alignments and competitive HiFi remapping. '
              'It does not constitute a complete repeat annotation or haplotype phasing analysis. '
              'Same-chromosome controls and overlapping scan windows are descriptive comparisons, not independent statistical replicates. '
              'Contact maps and underlying compressed count matrices are retained for further review.','']
    (a.out/'README.md').write_text('\n'.join(lines),encoding='utf-8')
    write(a.out/'sites.tsv',summaries,list(summaries[0]))
    write(a.out/'manual-cut-worksheet.tsv',actions,list(actions[0]))
    save(a.out/'provenance.json',manifest)


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    prep=sub.add_parser('prepare')
    for key in ('source','focus','final-dir','out'):prep.add_argument('--'+key,type=Path,required=True)
    branch=sub.add_parser('branch');branch.add_argument('--kind',choices=('labels','hifi','hic'),required=True)
    branch.add_argument('--job',type=Path,required=True);branch.add_argument('--out',type=Path,required=True)
    branch.add_argument('--threads',type=int,default=8);branch.add_argument('--masks',type=Path)
    rep=sub.add_parser('report');rep.add_argument('--manifest',type=Path,required=True)
    rep.add_argument('--packet',dest='packets',action='append',type=Path,required=True);rep.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.command=='prepare':prepare(a)
    elif a.command=='report':report(a)
    else:
        a.out.mkdir(exist_ok=False);job=json.loads(a.job.read_text())
        if a.kind=='hic':hic(job,a.out,a.threads,a.masks)
        elif a.kind=='hifi':hifi(job,a.out,a.threads,a.masks)
        else:labels(job,a.out,a.threads)
        save(a.out/'provenance.json',dict(assembly=job['assembly'],assessment_sha256=job['assessment_sha256'],
             coordinate_stage=job['coordinate_stage'],job_manifest_sha256=sha(a.job),
             script_sha256=sha(__file__),commands=COMMANDS,
             tools={name:subprocess.check_output([name,'--version'],text=True,errors='replace').splitlines()[0]
                    for name in ('samtools','minimap2')}))


if __name__=='__main__':main()
