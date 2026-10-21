"""General chromosome-junction assays shared by production and focused investigations.

Coordinates are zero-based, half-open. Unique means unique among retained alignments.
"""
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import re
import statistics
import subprocess
import sys
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

def refine_advice(event):
    """Positive placement evidence modifies review advice, never selects a cut.

    Keep this supplement independent of the older mapping/assessment helper so a
    review policy update does not invalidate whole-assembly HiFi mapping tasks.
    """
    reason='Multiple credible HiFi molecules reach unique sequence from both chromosome arms at this assay site; retain pending review of conflicting evidence.'
    supported=set()
    for target in [event]+event.get('alternatives',[])+event.get('cut_options',[]):
        evidence=target.get('junction_support',{})
        if evidence.get('credible_bridges') and target in event.get('cut_options',[]):
            target['automatic_cut_eligible']=False
        if (evidence.get('unique_arm_bridges') or 0)<2:continue
        supported.add(target['id'])
        target.update(site_verdict='REVIEW',site_reason=reason)
        if target.get('kind')=='regional_read_loss':
            target.update(rank_score=target.get('rank_score',0)-20,assessment='Review: credible molecules support this site',reason=reason)
    alternatives=event.get('alternatives',[])
    alternatives.sort(key=lambda t:(-t.get('rank_score',0),t.get('distance_bp',0),t['lo']))
    for rank,target in enumerate(alternatives,1):target['rank']=rank
    cut=event.get('recommended_cut',{})
    preferred=event.get('preferred_region',{})
    if cut.get('id') in supported or preferred.get('id') in supported:
        event.pop('recommended_cut',None);event.pop('preferred_region',None)
        event.update(recommendation='REVIEW',reason=reason)
    elif event['id'] in supported and not cut and not preferred:
        event.update(recommendation='REVIEW',reason=reason)


def credible_bridge(row):
    """Require confident placement and substantial high-identity read alignment."""
    aligned=int(row['aligned']);clipped=int(row['clipped'])
    return (int(row['mapq'])>=30 and int(row['same_haplotype_alternatives'])==0
            and float(row['identity'])>=.95 and aligned/max(1,aligned+clipped)>=.8)


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

def distance_matched_contacts(base, ls, rs, unique_bins, minimum_controls=5, cache=None):
    """Descriptive baseline from same-arm bin pairs at comparable distance/depth.

    Match genomic separation (10% or one bin) and per-bin diagonal contact depth
    (within a factor of two). No significance or independent-replicate claim.
    Every candidate bin pair must have controls; missing matches never become zero.
    """
    import numpy as np
    n=len(base);diagonal=np.diag(base).astype(float)
    cache={} if cache is None else cache
    pools=cache.get('pools',{})
    for chrom,bins in ([] if pools else unique_bins.items()):
        bins=sorted(i for i in bins if diagonal[i]>=5)
        for x,i in enumerate(bins):
            for j in bins[x+1:]:
                pools.setdefault(j-i,[]).append((i,j))
    cache['pools']=pools
    expected=0.;control_pairs=set();matched=0;total=len(ls)*len(rs)
    for i in ls:
        for j in rs:
            key=(i,j)
            if key in cache:
                value,controls=cache[key]
                if value is not None:expected+=value;matched+=1
                control_pairs.update(controls);continue
            distance=j-i;tolerance=max(1,math.ceil(distance*.1));values=[];local_controls=set()
            if min(diagonal[i],diagonal[j])<5:continue
            for d in range(max(1,distance-tolerance),distance+tolerance+1):
                for a,b in pools.get(d,[]):
                    if (a,b)==(i,j):continue
                    # Match each end separately; never use a chromosome-crossing pair.
                    if not (.5<=diagonal[a]/diagonal[i]<=2 and .5<=diagonal[b]/diagonal[j]<=2):continue
                    scale=math.sqrt(diagonal[i]*diagonal[j]/(diagonal[a]*diagonal[b]))
                    values.append(float(base[a,b])*scale);local_controls.add((a,b))
            value=sum(values)/len(values) if len(values)>=minimum_controls else None
            cache[key]=(value,local_controls)
            control_pairs.update(local_controls)
            if value is not None:expected+=value;matched+=1
    calibrated=bool(total and matched==total and expected>0)
    observed=int(base[np.ix_(ls,rs)].sum()) if total else None
    return dict(control_method='distance_and_diagonal_depth',matched_controls=len(control_pairs),
                matched_candidate_pairs=matched,candidate_pairs=total,
                calibration_status='measured' if calibrated else 'insufficient_matched_controls',
                expected_cross=expected if calibrated else None,
                observed_expected=observed/expected if calibrated else None)


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
        for lib in libs+(['POOLED'] if job.get('include_pooled',True) else []):
            base=sum((matrices[e['event'],l].astype(np.uint64) for l in libs),np.zeros((n,n),dtype=np.uint64)) if lib=='POOLED' else matrices[e['event'],lib]
            calibration_cache={}
            for res in job.get('resolutions',(10000,25000,50000,100000)):
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
                left=[i for i in sorted(unique_bins[e['left']]) if lo+(i+1)*resolution<=s['point'] and base[i,i]>=5]
                right=[i for i in sorted(unique_bins[e['right']]) if lo+i*resolution>=s['point'] and base[i,i]>=5]
                for width in (50000,100000,300000):
                    need=width//resolution;ls=left[-need:];rs=right[:need]
                    row=dict(event=e['event'],site=s['id'],point=s['point'],library=lib,anchor_bp=width,
                             left_bins=len(ls),right_bins=len(rs),status='insufficient_observable_unique_arm_sequence')
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
                        row.update(distance_matched_contacts(base,ls,rs,unique_bins,cache=calibration_cache))
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
    expected_replay=sum(len(rows) for e in old['events'] for rows in e.get('contact_profiles',{}).values())
    if expected_replay and not replayed:raise ValueError('No original contact measurements could be replayed')
    write(out/'contact_profile.tsv',profiles,['event','library','point','flank_bp','cross','left_within','right_within','ratio'])
    write(out/'unique_arm_contacts.tsv',anchors,['event','site','point','library','anchor_bp','status','left_bins','right_bins','left_lo','left_hi','right_lo','right_hi','left_intervals','right_intervals','cross','left_within','right_within','ratio','control_method','matched_controls','matched_candidate_pairs','candidate_pairs','calibration_status','expected_cross','observed_expected'])
    write(out/'controls.tsv',controls,['event','site','kind','library','point','flank_bp','chrom','cross','left_within','right_within','ratio'])
    write(out/'marginals.tsv',marginals,['event','library','lo','hi','local_contact_ends'])
    save(out/'summary.json',dict(branch='hic',assembly=job['assembly'],audit=dict(audit),maps=map_records,replayed_measurements=replayed,
         limitations=['Local marginals count contacts retained within the displayed window, not genome-wide coverage.',
         'Pooled measurements supplement separate libraries and are not an independent replicate.',
         'Overlapping scan positions and control windows are correlated; no independent-bin significance is claimed.']))
