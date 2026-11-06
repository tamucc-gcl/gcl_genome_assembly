"""Validate PSMC iterations, reproducibly resample blocks and export explicit curves."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import re


def read_fasta(path):
    name=None;seq=[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                if name is not None:yield name,''.join(seq)
                name=line[1:].split()[0];seq=[]
            elif line:seq.append(line)
    if name is not None:yield name,''.join(seq)


def resample(source,destination,sample,replicate,seed):
    blocks=list(read_fasta(source))
    if not blocks or any(not s or set(s)-set('KTN') for _,s in blocks):
        raise ValueError('Empty/invalid split PSMCFA')
    derived=int.from_bytes(hashlib.sha256(f'{seed}:{sample}:{replicate}'.encode()).digest()[:8],'big')
    rng=random.Random(derived);selection=[];target=sum(len(seq) for _,seq in blocks);total=0
    # Match the upstream length-targeted stopping rule, rather than assuming all
    # splitfa blocks have equal length. Only the RNG seed is supplied externally.
    while True:
        chosen=rng.randrange(len(blocks));length=len(blocks[chosen][1])
        deficit=target-total;overshoot=total+length-target
        if overshoot<=0 or (deficit>0 and overshoot<deficit):
            selection.append(chosen);total+=length
        if deficit>=0 and overshoot>=0:break
    with open(destination,'w') as f:
        for i,index in enumerate(selection):
            f.write(f'>bootstrap_{i+1}\n{blocks[index][1]}\n')
    result=dict(sample=sample,replicate=replicate,base_seed=seed,derived_seed=derived,
                rng='Python random.Random / MT19937',source_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),
                block_count=len(blocks),selected_block_indices=selection,target_bins=target,resampled_bins=total,
                method='Seeded splitfa block resampling with upstream nearest-total-length stopping rule; fit without native -b')
    Path('resampling.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def pattern_groups(pattern):
    if not re.fullmatch(r'\d+(?:\*\d+)?(?:\+\d+(?:\*\d+)?)*',pattern):raise ValueError('Invalid PSMC interval pattern')
    result=[]
    for token in pattern.split('+'):
        fields=[int(x) for x in token.split('*')]
        if any(x<=0 for x in fields):raise ValueError('Interval pattern values must be positive')
        result.extend([fields[-1]]*(fields[0] if len(fields)==2 else 1))
    if sum(result)>255:raise ValueError('PSMC supports at most 255 atomic intervals')
    return result


def parse(path):
    blocks=[];current=None;skip=None;pattern=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('MM'):
            found=re.search(r'skip:(\d+)',line)
            if found:skip=int(found[1])
            found=re.search(r'pattern:([^,\s]+)',line)
            if found:pattern=found[1]
        fields=line.split()
        if not fields:continue
        if fields[0]=='RD':
            if current is not None:raise ValueError('Unclosed PSMC iteration')
            current=dict(iteration=int(fields[1]),rows=[])
        elif fields[0] in ['TR','RS','DT'] and current is not None:
            values=[float(v) for v in fields[1:]]
            if any(not math.isfinite(v) for v in values):raise ValueError('Non-finite PSMC estimate')
            if fields[0]=='TR':current['theta']=values[0];current['rho']=values[1]
            elif fields[0]=='DT':current['dt']=values[0]
            else:current['rows'].append(values)
        elif fields[0]=='//' and current is not None:
            if not current.get('theta',0)>0 or not current['rows']:raise ValueError('Incomplete/invalid PSMC iteration')
            indices=[int(r[0]) for r in current['rows']]
            if indices!=list(range(len(indices))) or any(r[1]<0 or r[2]<=0 for r in current['rows']):
                raise ValueError('Invalid PSMC interval order or estimates')
            if any(b[1]<a[1] for a,b in zip(current['rows'],current['rows'][1:])):
                raise ValueError('Non-monotone PSMC times')
            blocks.append(current);current=None
    if current is not None:raise ValueError('Truncated PSMC output: missing closing //')
    if not blocks or not skip or not pattern:raise ValueError('PSMC output lacks fitted iterations/metadata')
    if len({b['iteration'] for b in blocks})!=len(blocks):raise ValueError('Duplicate PSMC iteration')
    groups=pattern_groups(pattern)
    if any(len(b['rows'])!=sum(groups) for b in blocks):raise ValueError('PSMC interval count disagrees with pattern')
    return dict(skip=skip,pattern=pattern,blocks=blocks)


def export(files,out,iteration,bin_size,mutation_rate=None,generation_time=None,show_last=False,expected_bootstraps=0):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    names=[Path(f).name for f in files]
    expected={'observed.psmc'} | {f'bootstrap_{i:04d}.psmc' for i in range(1,expected_bootstraps+1)}
    if len(names)!=len(set(names)) or set(names)!=expected:raise ValueError('Missing/unexpected observed/bootstrap fit set')
    if bin_size<=0:raise ValueError('Bin size must be positive')
    if (mutation_rate is None)!=(generation_time is None):raise ValueError('Both scaling values are required')
    if mutation_rate is not None and (not math.isfinite(mutation_rate) or not math.isfinite(generation_time) or min(mutation_rate,generation_time)<=0):
        raise ValueError('Scaling must be finite and positive')
    rows=[];fits=[];interval_qc=[]
    for path in sorted(files,key=lambda f:Path(f).name!='observed.psmc'):
        data=parse(path);chosen=[b for b in data['blocks'] if b['iteration']==iteration]
        if len(chosen)!=1:raise ValueError('Requested plotting iteration absent: '+str(path))
        block=chosen[0];theta=block['theta']/(bin_size*data['skip'])
        groups=pattern_groups(data['pattern']);last_start=sum(groups[:-1]);start=0
        for group,span in enumerate(groups):
            recombinations=sum(r[3] for r in block['rows'][start:start+span])
            interval_qc.append(dict(fit=Path(path).stem,parameter_group=group,first_atomic_interval=start,
                last_atomic_interval=start+span-1,expected_recombinations=recombinations,
                below_upstream_advisory_ten='yes' if recombinations<10 else 'no'))
            start+=span
        for r in block['rows']:
            index=int(r[0])
            if not show_last and index>=last_start:continue
            raw_t=r[1]+block.get('dt',0);lam=r[2]
            if mutation_rate is None:x=raw_t*theta;y=lam*theta
            else:
                n0=theta/(4*mutation_rate);x=2*n0*raw_t*generation_time;y=n0*lam
            rows.append(dict(fit=Path(path).stem,kind='observed' if Path(path).name=='observed.psmc' else 'bootstrap',
                             iteration=iteration,interval=index,time=x,population=y,raw_time=raw_t,relative_size=lam,
                             theta_per_base=theta,expected_recombinations=r[3],
                             scaling='unscaled' if mutation_rate is None else 'years_Ne'))
        fits.append(dict(file=Path(path).name,sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                         iteration=iteration,last_complete_iteration=max(b['iteration'] for b in data['blocks']),
                         pattern=data['pattern'],internal_skip=data['skip']))
    with (out/'curves.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    with (out/'interval_qc.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(interval_qc[0]),delimiter='\t');w.writeheader();w.writerows(interval_qc)
    (out/'fit_manifest.json').write_text(json.dumps(dict(fits=fits,bin_size=bin_size,mutation_rate=mutation_rate,
        generation_time=generation_time,show_last_interval=show_last,
        method='PSMC upstream scaling; unscaled time=t*theta/bin and size=lambda*theta/bin; terminal free-parameter block omitted by default'),indent=2)+'\n')
    return rows



def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('resample');q.add_argument('--input',required=True);q.add_argument('--out',required=True);q.add_argument('--sample',required=True)
    q.add_argument('--replicate',type=int,required=True);q.add_argument('--seed',type=int,required=True)
    q=sub.add_parser('validate');q.add_argument('--input',required=True);q.add_argument('--last-iteration',type=int,required=True)
    q=sub.add_parser('export');q.add_argument('--files',nargs='+',required=True);q.add_argument('--out',required=True)
    q.add_argument('--iteration',type=int,required=True);q.add_argument('--bin-size',type=int,required=True)
    q.add_argument('--mutation-rate',type=float);q.add_argument('--generation-time',type=float)
    q.add_argument('--show-last',action='store_true');q.add_argument('--expected-bootstraps',type=int,default=0)
    a=p.parse_args()
    if a.command=='resample':resample(a.input,a.out,a.sample,a.replicate,a.seed)
    elif a.command=='validate':
        data=parse(a.input)
        if max(b['iteration'] for b in data['blocks'])<a.last_iteration:raise ValueError('PSMC stopped before requested iteration')
    elif a.command=='export':export(a.files,a.out,a.iteration,a.bin_size,a.mutation_rate,a.generation_time,a.show_last,a.expected_bootstraps)


if __name__=='__main__':main()
