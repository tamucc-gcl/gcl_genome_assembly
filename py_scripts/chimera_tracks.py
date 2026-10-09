"""Diagnostic chromosome tracks; assignments never authorize cuts."""
import re
import gzip
from collections import defaultdict,Counter
from bisect import bisect_left,bisect_right


class LabelIndex:
    """Disjoint chromosome runs and logarithmic coordinate/edge lookup.

    Sweep overlapping input labels once. Conflicting chromosomes remain unassigned;
    adjacent equivalent assignments merge without bridging unlabelled sequence.
    """
    def __init__(self, labels):
        self.plain={};self.runs={};self.starts={};self.edges={};self.raw_edges={}
        for target,value in labels.items():
            if isinstance(value,str):
                self.plain[target]=value;continue
            events=defaultdict(Counter)
            for b in value or []:
                if b['hi']<=b['lo']:continue
                events[b['lo']][b['chrom']]+=1;events[b['hi']][b['chrom']]-=1
            coordinates=sorted(events);active=Counter();runs=[]
            for i,lo in enumerate(coordinates[:-1]):
                active.update(events[lo]);active=Counter({c:n for c,n in active.items() if n>0})
                chrom=next(iter(active)) if len(active)==1 else None
                hi=coordinates[i+1]
                if chrom is not None:
                    if runs and runs[-1][1]==lo and runs[-1][2]==chrom:runs[-1]=(runs[-1][0],hi,chrom)
                    else:runs.append((lo,hi,chrom))
            self.runs[target]=runs;self.starts[target]=[r[0] for r in runs]
            self.edges[target]=sorted({edge for lo,hi,_ in runs for edge in (lo,hi)})
            self.raw_edges[target]=coordinates

    def at(self,target,coordinate):
        if coordinate is None:return None
        if target in self.plain:return self.plain[target]
        i=bisect_right(self.starts.get(target,[]),coordinate)-1
        runs=self.runs.get(target,[])
        return runs[i][2] if i>=0 and coordinate<runs[i][1] else None

    def inside(self,target,lo,hi,compact):
        edges=(self.edges if compact else self.raw_edges).get(target,[])
        return edges[bisect_right(edges,lo):bisect_left(edges,hi)]


def chromosome_at(labels, target, coordinate):
    if coordinate is None:return None
    value=labels.get(target)
    if isinstance(value,str):return value
    matches={b['chrom'] for b in (value or []) if b['lo']<=coordinate<b['hi']}
    return next(iter(matches)) if len(matches)==1 else None


def augment_peer_labels(peer):
    """Label blocks inside composites from their own reference PAF, never whole composites."""
    path=peer.get('reference_paf','')
    if not path or path.split('/')[-1]=='NO_PAF':return
    from chimera_intervals import regions
    hits=defaultdict(list)
    opener=gzip.open if path.endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8') as handle:
        for line in handle:
            f=line.rstrip().split('\t')
            if len(f)<12:raise ValueError('Malformed peer reference PAF')
            chrom=peer.get('reference_labels',{}).get(f[5].split('#')[-1])
            if chrom and f[0] not in peer['chromosome_labels'] and int(f[11])>=20 and int(f[9])/max(1,int(f[10]))>=.9:
                hits[f[0]].append((int(f[2]),int(f[3]),chrom))
    for scaffold,values in hits.items():
        # Require 100 kb of union evidence per chromosome, not one uninterrupted
        # 100 kb interval. Ambiguity and unaligned gaps remain explicitly unlabelled.
        pieces=regions(values,1);support=Counter()
        for lo,hi,chrom in pieces:support[chrom]+=hi-lo
        peer['chromosome_labels'][scaffold]=[dict(lo=lo,hi=hi,chrom=chrom) for lo,hi,chrom in pieces if support[chrom]>=100000]
    peer['block_label_source']='Own reference PAF; >=100 kb union support per chromosome; ambiguous/unmapped bases unlabelled; MAPQ >=20, identity >=90%'

def blocks(fields,compact=False):
    cigar=next((x[5:] for x in fields[12:] if x.startswith('cg:Z:')),None)
    if not cigar:return []
    ops=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if ''.join(n+op for n,op in ops)!=cigar:raise ValueError('Malformed CIGAR')
    q=int(fields[2]) if fields[4]=='+' else int(fields[3]);target=int(fields[7]);result=[]
    for n,op in ops:
        n=int(n)
        if op in 'M=X':
            lo,hi=(q,q+n) if fields[4]=='+' else (q-n,q)
            if compact and result and result[-1]['t']+result[-1]['hi']-result[-1]['lo']==target and (
                    result[-1]['hi']==lo if fields[4]=='+' else result[-1]['lo']==hi):
                result[-1]['lo']=min(result[-1]['lo'],lo);result[-1]['hi']=max(result[-1]['hi'],hi)
            else:
                result.append(dict(lo=lo,hi=hi,t=target,strand=fields[4],target=fields[5],length=int(fields[6]),
                    mapq=int(fields[11]),identity=int(fields[9])/max(1,int(fields[10]))))
        if op in 'MI=X':q+=n if fields[4]=='+' else -n
        if op in 'MDN=X':target+=n
    if q!=(int(fields[3]) if fields[4]=='+' else int(fields[2])) or target!=int(fields[8]):raise ValueError('CIGAR endpoint mismatch')
    return result


def position(block,query):
    return block['t']+(query-block['lo'] if block['strand']=='+' else block['hi']-1-query)


def tracks(hits,labels,mode="locus_unique",min_mapq=30,minimum_target=5000000,compact=False):
    label_index=LabelIndex(labels)
    segments=[b for f in hits for b in blocks(f,compact)]
    events=defaultdict(list)
    for i,b in enumerate(segments):events[b['lo']].append((i,1));events[b['hi']].append((i,-1))
    for b in segments:
        for target_edge in label_index.inside(b['target'],b['t'],b['t']+b['hi']-b['lo'],compact):
            q=b['lo']+target_edge-b['t'] if b['strand']=='+' else b['hi']-(target_edge-b['t'])
            events[q]
    active=set();result=[];coordinates=sorted(events)
    for index,lo in enumerate(coordinates[:-1]):
        for i,change in events[lo]:
            if change==1:active.add(i)
            else:active.discard(i)
        hi=coordinates[index+1]
        qualified=[segments[i] for i in active if segments[i]['mapq']>=min_mapq and segments[i]['identity']>=.9 and
                   segments[i]['length']>=minimum_target and segments[i]['target'] in labels]
        chrom=None;status='unaligned_or_filtered'
        if qualified:
            best=max(qualified,key=lambda b:(b['mapq'],b['identity']))
            mid=(lo+hi)//2
            best_chrom=label_index.at(best['target'],position(best,mid))
            competing=[segments[i] for i in active if segments[i]['identity']>=.95*best['identity'] and
                (segments[i]['target']!=best['target'] or segments[i]['strand']!=best['strand'] or abs(position(segments[i],mid)-position(best,mid))>1000)]
            if competing and (mode=='locus_unique' or any(label_index.at(b['target'],position(b,mid))!=best_chrom for b in competing)):
                status='ambiguous'
            else:
                chrom=best_chrom;status=('chromosome_assigned_locus_ambiguous' if competing else 'assigned') if chrom else 'unlabelled_target_block'
        if compact and result and result[-1]['hi']==lo and result[-1]['chrom']==chrom and result[-1]['status']==status:
            result[-1]['hi']=hi
        else:result.append(dict(lo=lo,hi=hi,chrom=chrom,status=status))
    return result


def bins(track,length,width=50000):
    result=[]
    for lo in range(0,length,width):
        hi=min(length,lo+width);counts=Counter();ambiguous=0
        for b in track:
            size=max(0,min(hi,b['hi'])-max(lo,b['lo']))
            if b['chrom']:counts[b['chrom']]+=size
            elif b['status']=='ambiguous':ambiguous+=size
        total=sum(counts.values());chrom,dominant=max(counts.items(),key=lambda x:x[1]) if counts else (None,0)
        result.append(dict(lo=lo,hi=hi,aligned_bp=total,ambiguous_bp=ambiguous,dominance=dominant/max(1,total),
            chrom=chrom if total>=10000 and dominant/max(1,total)>=.9 else None,counts=dict(counts)))
    return result


def side_summary(track,tiles,lo,hi):
    counts=Counter()
    for b in track:
        if b['chrom']:counts[b['chrom']]+=max(0,min(hi,b['hi'])-max(lo,b['lo']))
    total=sum(counts.values());chrom,dominant=max(counts.items(),key=lambda x:x[1]) if counts else (None,0)
    informative=[b for b in tiles if b['lo']>=lo and b['hi']<=hi and chrom is not None and b['chrom']==chrom]
    return dict(chrom=chrom,aligned_bp=total,coverage=total/max(1,hi-lo),dominance=dominant/max(1,total),
        informative_bins=len(informative),qualified=total>=100000 and dominant/max(1,total)>=.9 and len(informative)>=3)

