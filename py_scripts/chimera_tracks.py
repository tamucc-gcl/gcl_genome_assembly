"""Diagnostic chromosome tracks; assignments never authorize cuts."""
import re
from collections import defaultdict,Counter

def blocks(fields):
    cigar=next((x[5:] for x in fields[12:] if x.startswith('cg:Z:')),None)
    if not cigar:return []
    ops=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if ''.join(n+op for n,op in ops)!=cigar:raise ValueError('Malformed CIGAR')
    q=int(fields[2]) if fields[4]=='+' else int(fields[3]);target=int(fields[7]);result=[]
    for n,op in ops:
        n=int(n)
        if op in 'M=X':
            lo,hi=(q,q+n) if fields[4]=='+' else (q-n,q)
            result.append(dict(lo=lo,hi=hi,t=target,strand=fields[4],target=fields[5],length=int(fields[6]),
                mapq=int(fields[11]),identity=int(fields[9])/max(1,int(fields[10]))))
        if op in 'MI=X':q+=n if fields[4]=='+' else -n
        if op in 'MDN=X':target+=n
    if q!=(int(fields[3]) if fields[4]=='+' else int(fields[2])) or target!=int(fields[8]):raise ValueError('CIGAR endpoint mismatch')
    return result


def position(block,query):
    return block['t']+(query-block['lo'] if block['strand']=='+' else block['hi']-1-query)


def tracks(hits,labels,mode="locus_unique"):
    segments=[b for f in hits for b in blocks(f)]
    events=defaultdict(list)
    for i,b in enumerate(segments):events[b['lo']].append((i,1));events[b['hi']].append((i,-1))
    active=set();result=[];coordinates=sorted(events)
    for index,lo in enumerate(coordinates[:-1]):
        for i,change in events[lo]:
            if change==1:active.add(i)
            else:active.discard(i)
        hi=coordinates[index+1]
        qualified=[segments[i] for i in active if segments[i]['mapq']>=30 and segments[i]['identity']>=.9 and
                   segments[i]['length']>=5000000 and segments[i]['target'] in labels]
        chrom=None;status='unaligned_or_filtered'
        if qualified:
            best=max(qualified,key=lambda b:(b['mapq'],b['identity']))
            mid=(lo+hi)//2
            competing=[segments[i] for i in active if segments[i]['identity']>=.95*best['identity'] and
                (segments[i]['target']!=best['target'] or segments[i]['strand']!=best['strand'] or abs(position(segments[i],mid)-position(best,mid))>1000)]
            if competing and (mode=='locus_unique' or any(labels.get(b['target'])!=labels[best['target']] for b in competing)):
                status='ambiguous'
            else:
                chrom=labels[best['target']];status='chromosome_assigned_locus_ambiguous' if competing else 'assigned'
        result.append(dict(lo=lo,hi=hi,chrom=chrom,status=status))
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

