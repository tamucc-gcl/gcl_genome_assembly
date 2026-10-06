"""Measured junction/control assays. Zero counts never imply an informative assay."""
from collections import defaultdict, Counter
import csv
import gzip
import math
import re
import statistics
import hashlib
import bisect


def agp(path):
    gaps, placements, lengths = [], defaultdict(list), {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip().split('\t')
            if len(f) != 9:
                raise ValueError('Invalid AGP row')
            lo, hi = int(f[1])-1, int(f[2])
            lengths[f[0]] = max(lengths.get(f[0], 0), hi)
            if f[4] in ('N', 'U'):
                if hi-lo != int(f[5]):
                    raise ValueError('AGP gap length mismatch')
                gaps.append(dict(scaffold=f[0], lo=lo, hi=hi))
            else:
                start, end = int(f[6])-1, int(f[7])
                if hi-lo != end-start or f[8] not in ('+', '-'):
                    raise ValueError('AGP component is not exactly liftable')
                placements[f[5]].append((f[0], lo, hi, start, end, f[8]))
    return gaps, placements, lengths


def select_controls(intervals, gaps, lengths, count=12, flank=250000):
    """Match gap length and complete flank geometry; avoid every suspect transition."""
    controls, links = {}, {}
    for key, interval in list(intervals.items()):
        compatible = [g for g in gaps if g['scaffold']==interval['scaffold'] and
                      interval['lo'] <= g['lo'] < g['hi'] <= interval['hi']]
        if len(compatible) != 1:
            links[key] = []; continue
        focal = compatible[0]
        interval['gap'] = focal
        size = focal['hi']-focal['lo']
        available = [g for g in gaps if .5*size <= g['hi']-g['lo'] <= 2*size and
                     g['lo'] >= flank and lengths[g['scaffold']]-g['hi'] >= flank and
                     not any(g['scaffold']==s['scaffold'] and g['lo'] < s['hi']+flank and
                             g['hi'] > s['lo']-flank for s in intervals.values())]
        available.sort(key=lambda g:(g['scaffold'] != focal['scaffold'],
                                     abs((g['lo']+g['hi'])-(focal['lo']+focal['hi'])), g['scaffold'], g['lo']))
        chosen = []
        for g in available:
            if any(g['scaffold']==x['scaffold'] and abs(g['lo']-x['lo']) < 2*flank for x in chosen):
                continue
            chosen.append(g)
            if len(chosen)==count:
                break
        links[key] = []
        for g in chosen:
            cid = 'control_'+hashlib.sha256(('%s:%d:%d' % (g['scaffold'],g['lo'],g['hi'])).encode()).hexdigest()[:20]
            controls[cid] = dict(scaffold=g['scaffold'], lo=g['lo'], hi=g['hi'], gap=g, role='control')
            links[key].append(cid)
    return controls, links


def nearby_hypotheses(intervals, gaps, radius=250000, maximum=8):
    """Discover review-only gap hypotheses; do not replace a transition by its nearest gap."""
    result={}
    localized=set()
    for key,interval in intervals.items():
        compatible=[g for g in gaps if g['scaffold']==interval['scaffold'] and interval['lo']<=g['lo']<g['hi']<=interval['hi']]
        if len(compatible)==1:localized.add((compatible[0]['scaffold'],compatible[0]['lo'],compatible[0]['hi']))
    for key,interval in intervals.items():
        choices=[g for g in gaps if g['scaffold']==interval['scaffold'] and g['lo']<interval['hi']+radius and g['hi']>interval['lo']-radius and
                 (g['scaffold'],g['lo'],g['hi']) not in localized]
        choices.sort(key=lambda g:(min(abs(g['lo']-interval['hi']),abs(g['hi']-interval['lo'])),g['lo']))
        for g in choices[:maximum]:
            hid='hypothesis_'+hashlib.sha256(('%s:%d:%d'%(g['scaffold'],g['lo'],g['hi'])).encode()).hexdigest()[:20]
            value=result.setdefault(hid,dict(scaffold=g['scaffold'],lo=g['lo'],hi=g['hi'],gap=g,role='review_hypothesis',source_candidates=[]))
            value['source_candidates'].append(key)
    return result


def qualified_molecule(fields):
    if len(fields)<11 or int(fields[1])&(4|256|512|1024|2048) or int(fields[4])<30:
        return False
    operations=re.findall(r'(\d+)([MIDNSHP=X])',fields[5])
    if ''.join(n+op for n,op in operations)!=fields[5]:raise ValueError('Malformed SAM CIGAR')
    aligned=sum(int(n) for n,op in operations if op in 'MI=X')
    nm=next((int(tag[5:]) for tag in fields[11:] if tag.startswith('NM:i:')),None)
    if nm is None:
        if any(op=='M' for n,op in operations):return False
        nm=sum(int(n) for n,op in operations if op in 'IDX')
    return aligned>0 and nm>=0 and nm/aligned<=.1


def sam_measure(path, lo, hi, anchor=1000):
    """Collapse molecules; reject large deletions/skips as uninterrupted spans."""
    left, right, spanning = set(), set(), set()
    left_depth, right_depth = [0]*anchor, [0]*anchor
    seen = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            f = line.rstrip().split('\t')
            if not qualified_molecule(f) or f[0] in seen:
                continue
            seen.add(f[0])
            pos, blocks, bad = int(f[3])-1, [], False
            cigar = re.findall(r'(\d+)([MIDNSHP=X])', f[5])
            if ''.join(n+op for n,op in cigar) != f[5]:
                raise ValueError('Malformed SAM CIGAR')
            for n, op in cigar:
                n=int(n)
                if op in 'M=X':
                    blocks.append((pos,pos+n)); pos+=n
                elif op in 'DN':
                    if n>50 and pos<hi+anchor and pos+n>lo-anchor:
                        bad=True
                    pos+=n
            coverage=[]
            for start, end, depth in ((lo-anchor,lo,left_depth),(hi,hi+anchor,right_depth)):
                bases=set()
                for begin, finish in blocks:
                    bases.update(range(max(0,begin-start), min(anchor,finish-start)))
                coverage.append(len(bases))
                for base in bases: depth[base]+=1
            if coverage[0]>=anchor: left.add(f[0])
            if coverage[1]>=anchor: right.add(f[0])
            if not bad and min(coverage)>=anchor: spanning.add(f[0])
    return dict(left_molecules=len(left), right_molecules=len(right), spanning=len(spanning),
                left_median_depth=statistics.median(left_depth), right_median_depth=statistics.median(right_depth),
                left_covered_fraction=sum(x>=5 for x in left_depth)/anchor,
                right_covered_fraction=sum(x>=5 for x in right_depth)/anchor)


def project_query(f,query):
    cigar=next((x[5:] for x in f[12:] if x.startswith('cg:Z:')),None)
    if not cigar:return None
    q=int(f[2]) if f[4]=='+' else int(f[3]);t=int(f[7])
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cigar):
        n=int(n)
        if op in 'M=X':
            start,end=(q,q+n) if f[4]=='+' else (q-n,q)
            if start<=query<end:return t+(query-q if f[4]=='+' else q-1-query)
        if op in 'MI=X':q+=n if f[4]=='+' else -n
        if op in 'MDN=X':t+=n
    return None


def flank_placement(hits, lo, hi, minimum_target=5000000):
    qualified=[]
    for f in hits:
        if int(f[11])<30 or int(f[6])<minimum_target or int(f[9])/max(1,int(f[10]))<.9:
            continue
        covered=max(0,min(hi,int(f[3]))-max(lo,int(f[2])))
        if covered < .9*(hi-lo):
            continue
        tags={x[:2]:x[5:] for x in f[12:]}
        cigar=tags.get('cg')
        if not cigar:
            continue
        operations=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
        if ''.join(n+op for n,op in operations)!=cigar:
            raise ValueError('Malformed PAF CIGAR')
        if sum(int(n) for n,op in operations if op in 'MI=X') != int(f[3])-int(f[2]) or sum(int(n) for n,op in operations if op in 'MDN=X') != int(f[8])-int(f[7]):
            raise ValueError('PAF CIGAR endpoints do not agree')
        qpos=int(f[2]) if f[4]=='+' else int(f[3])
        aligned=0
        for n,op in operations:
            n=int(n)
            if op in 'M=X':
                start,end=(qpos,qpos+n) if f[4]=='+' else (qpos-n,qpos)
                aligned+=max(0,min(hi,end)-max(lo,start))
            if op in 'MI=X': qpos+=n if f[4]=='+' else -n
        if aligned >= .9*(hi-lo): qualified.append(f)
    if not qualified:
        return None
    targets={(f[5],f[4]) for f in qualified}
    # Retain secondary/repeat alternatives even if their MAPQ is zero.
    best=max(qualified,key=lambda f:int(f[9]))
    midpoint=(lo+hi)//2
    best_position=project_query(best,midpoint)
    alternatives=[]
    for f in hits:
        if f is best:continue
        covered=max(0,min(hi,int(f[3]))-max(lo,int(f[2])))
        if covered<.9*(hi-lo) or int(f[9])/max(1,int(f[10])) < .95*int(best[9])/max(1,int(best[10])):continue
        position=project_query(f,midpoint)
        if (f[5],f[4])!=(best[5],best[4]) or position is None or best_position is None or abs(position-best_position)>1000:
            alternatives.append(f)
    return best if len(targets)==1 and not alternatives else None


def peer_relationship(hits, lo, hi, anchor=50000, minimum_target=5000000):
    left=flank_placement(hits,lo-anchor,lo,minimum_target)
    right=flank_placement(hits,hi,hi+anchor,minimum_target)
    if left is None or right is None:
        return 'uninformative', None, None
    if left[5]!=right[5]:
        return 'separate_scaffolds', left, right
    x,y=project_query(left,lo-1),project_query(right,hi)
    ordered=x is not None and y is not None and ((y>x) if left[4]=='+' else (y<x))
    if left[4]==right[4] and ordered and abs(y-x)<=hi-lo+10000:
        return 'continuous_context', left, right
    return 'uninformative', left, right


def scan_contacts(path, placements, intervals, libraries, flank=250000):
    """One pairs pass, preserving library identity and cross/within geometry."""
    windows=defaultdict(list)
    counts=defaultdict(Counter)
    prefixes=[r['qname_prefix'] for r in libraries]
    if not prefixes or any(not p for p in prefixes) or any(p.startswith(q) for i,p in enumerate(prefixes) for j,q in enumerate(prefixes) if i!=j):
        raise ValueError('Missing or overlapping Hi-C read-set prefixes')
    for key,r in intervals.items():
        for side,start,end in (('left',r['lo']-flank,r['lo']),('right',r['hi'],r['hi']+flank)):
            for b in range(max(0,start)//flank,(end-1)//flank+1):
                windows[r['scaffold'],b].append((key,side,start,end))
    def locate(chrom,pos):
        found=[v for v in placements.get(chrom,[]) if v[3]<=pos<v[4]]
        if len(found)!=1: return None
        v=found[0]
        return (v[0],v[1]+pos-v[3] if v[5]=='+' else v[2]-1-(pos-v[3]))
    total, audit = Counter(), Counter()
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip(): continue
            f=line.rstrip().split('\t'); audit['input_pairs']+=1
            if len(f)<8: raise ValueError('Malformed pairs row')
            if f[7]!='UU': continue
            ids={r['library_id'] for r in libraries if f[0].startswith(r['qname_prefix'])}
            if len(ids)!=1: raise ValueError('Unassigned or ambiguous Hi-C library')
            library=next(iter(ids)); total[library]+=1
            x,y=locate(f[1],int(f[2])-1),locate(f[3],int(f[4])-1)
            if x is None or y is None: audit['unliftable_pairs']+=1;continue
            ends=[]
            for chrom,pos in (x,y):
                ends.append({(key,side) for key,side,start,end in windows[chrom,pos//flank] if start<=pos<end})
            for key,side in ends[0]&ends[1]: counts[library,key][side+'_within']+=1
            for tags in ends:
                for key,side in tags: counts[library,key][side+'_ends']+=1
            crossed={key for key,side in ends[0] if (key,'right' if side=='left' else 'left') in ends[1]}
            for key in crossed: counts[library,key]['cross']+=1
    return counts, total, audit


def contact_ratio(counts):
    # Same geometry for each tested gap and length-matched controls.
    return counts.get('cross',0)/max(1,math.sqrt(counts.get('left_within',0)*counts.get('right_within',0)))


def informative_hifi(value):
    return min(value['left_molecules'],value['right_molecules'])>=10 and min(value['left_covered_fraction'],value['right_covered_fraction'])>=.95


def continuity_grid(path, lo, hi, anchor=1000, step=1000):
    """Narrow molecule spans across a bounded grid, not one broad-core bridge count."""
    points=sorted(set([lo,hi]+list(range(lo,hi+1,step))))
    counts=[0]*len(points);seen=set()
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):continue
            f=line.rstrip().split('\t')
            if not qualified_molecule(f) or f[0] in seen:continue
            seen.add(f[0]);pos=int(f[3])-1;indices=set()
            for n,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]):
                n=int(n)
                if op in 'M=X':
                    begin=bisect.bisect_left(points,pos+anchor)
                    end=bisect.bisect_right(points,pos+n-anchor)
                    indices.update(range(begin,end));pos+=n
                elif op in 'DN':pos+=n
            for i in indices:counts[i]+=1
    return dict(step_bp=step,anchor_bp=anchor,minimum_molecules=min(counts),
                supported_fraction=sum(n>=2 for n in counts)/len(counts),
                probes=[dict(cut_bp=point,molecules=n) for point,n in zip(points,counts)])


def compare_controls(candidate, controls, counts, libraries):
    """Explicit heuristic thresholds; an informative control is never assumed intact."""
    good=[c for c in controls if c.get('peer_continuous_individuals',0)>=3 and
          informative_hifi(c['hifi']) and c['hifi']['spanning']>=2 and
          all(.5<=c['hifi'][side+'_molecules']/max(1,candidate['hifi'][side+'_molecules'])<=2 for side in ('left','right'))]
    hifi_ok=informative_hifi(candidate['hifi'])
    results=[]
    for lib in sorted(libraries):
        focal=counts[lib,candidate['key']]
        matched=[c for c in good if min(counts[lib,c['key']].get('left_within',0),counts[lib,c['key']].get('right_within',0))>=100 and
                 .5 <= counts[lib,c['key']].get('left_within',0)/max(1,focal.get('left_within',0)) <=2 and
                 .5 <= counts[lib,c['key']].get('right_within',0)/max(1,focal.get('right_within',0)) <=2]
        informative=min(focal.get('left_within',0),focal.get('right_within',0))>=100 and len(matched)>=5
        ratios=[contact_ratio(counts[lib,c['key']]) for c in matched]
        minimum=min(ratios) if ratios else 0
        # Add a three-count allowance so zero observations cannot imply perfect certainty.
        upper=(focal.get('cross',0)+3)/max(1,math.sqrt(focal.get('left_within',0)*focal.get('right_within',0)))
        loss=informative and minimum>0 and upper < .1*minimum
        results.append(dict(library=lib, matched_controls=len(matched), informative=informative,
                            ratio=contact_ratio(focal), upper_count_allowance_ratio=upper,
                            minimum_control_ratio=minimum, support_loss=loss, raw_counts=dict(focal)))
    return dict(hifi_informative=hifi_ok, matched_controls_pass=hifi_ok and len(good)>=5,
                hic_informative=bool(results) and all(r['informative'] for r in results),
                hic_support_loss=bool(results) and all(r['support_loss'] for r in results), libraries=results)
