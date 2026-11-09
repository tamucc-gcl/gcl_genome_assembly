"""Removable read-support audit of one individual's actual masked PSMC input."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from psmc_reliability import fasta,read_support,write_table


def run(results,sample,out,sites=500):
    import numpy as np
    import pysam
    results=Path(results);out=Path(out)
    if out.exists():raise ValueError('Use a new read-support output directory')
    plan=json.loads((results/'post_assembly/psmc/psmc_plan.json').read_text())
    job=next(j for j in plan['jobs'] if j['sample']==sample)
    base=results/'post_assembly/psmc'/str(job['taxid'])/sample
    qc=json.loads((base/'input/input_qc.json').read_text())
    if sites<1:raise ValueError('Positive site count required')
    reference=job['reference'];mapping=results/'bam/hifi/final'
    bam=mapping/(reference['id']+'.sorted.bam')
    digest=mapping/(reference['id']+'.reference.sha256')
    if not bam.is_file() or not Path(str(bam)+'.bai').is_file():raise ValueError('Published final HiFi mapping/index missing: '+str(bam))
    if digest.read_text().split()[0]!=reference['sha256'] or qc['reference_sha256']!=reference['sha256']:
        raise ValueError('Final HiFi mapping and PSMC input references differ')
    rng=np.random.default_rng(1729);reservoir={};counts={};out.mkdir(parents=True)
    sha=hashlib.sha256()
    with open(reference['fasta'],'rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):sha.update(block)
    if sha.hexdigest()!=reference['sha256']:raise ValueError('Selected final reference content changed')
    with pysam.FastaFile(reference['fasta']) as fa:
        lengths=dict(zip(fa.references,fa.lengths))
        for name,seq in fasta(base/'input/diploid.fa.gz'):
            if len(seq)!=lengths[name]:raise ValueError('Consensus coordinate mismatch')
            for lo in range(0,len(seq),1000000):
                hi=min(len(seq),lo+1000000);calls=np.frombuffer(seq[lo:hi],dtype='S1')
                ref=np.frombuffer(fa.fetch(name,lo,hi).upper().encode(),dtype='S1')
                masks=dict(heterozygous=np.isin(calls,[c.encode() for c in 'MRWSYK']),
                    homozygous_reference=(calls==ref)&np.isin(calls,[b'A',b'C',b'G',b'T']))
                for group,mask in masks.items():
                    positions=np.flatnonzero(mask);counts[group]=counts.get(group,0)+len(positions);keys=rng.random(len(positions))
                    if len(keys)>sites:
                        pick=np.argpartition(keys,sites-1)[:sites];positions=positions[pick];keys=keys[pick]
                    values=reservoir.get(group,[])+[(float(k),name,int(lo+i),ref[i].decode(),calls[i].decode()) for k,i in zip(keys,positions)]
                    reservoir[group]=sorted(values)[:sites]
    if counts.get('heterozygous',0)!=qc['heterozygous_snps']:raise ValueError('Sampled input does not reproduce heterozygous-site count')
    selected=[dict(group=group,scaffold=n,pos_1based=p+1,reference=r,input=g)
        for group,rows in reservoir.items() for _,n,p,r,g in sorted(rows,key=lambda r:(r[1],r[2]))]
    evidence=list(read_support(bam,selected,lengths,('input',)))
    write_table(out/'sites.tsv',evidence,['group','scaffold','pos_1based','reference','input','molecules','A','C','G','T','input_support'])
    groups=[]
    for group in reservoir:
        rows=[r for r in evidence if r['group']==group];observable=[r for r in rows if r['molecules']>=10]
        groups.append(dict(group=group,population_sites=counts[group],sampled=len(rows),observable=len(observable),
            supported=sum(r['input_support']=='supported' for r in rows),insufficient=len(rows)-len(observable)))
    result=dict(sample=sample,mode=job['mode'],reference_sha256=reference['sha256'],bam=str(bam),seed=1729,groups=groups,
        criteria='MAPQ/baseQ >=20; >=10 distinct observed molecules; >=90% proposed allele support; each het allele >=20%; BAQ disabled')
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    md=[f'# PSMC input read support: {sample}','','Direct read support from the published, reference-verified final HiFi BAM.','',
        '| Site group | Population sites | Sampled | Observable | Supported | Insufficient reads |','| --- | ---: | ---: | ---: | ---: | ---: |']
    for g in groups:md.append(f"| {g['group']} | {g['population_sites']} | {g['sampled']} | {g['observable']} | {g['supported']} | {g['insufficient']} |")
    md+=['','[Per-site molecule counts](sites.tsv) · [Provenance](summary.json)','',result['criteria'],
        'These stratified samples assess read consistency, not genome-wide accuracy. Read mapping uses one haplotype and the reads contributed to assembly. No missing bases are filled, no genotypes changed, and no additional calls or fits run.']
    (out/'README.md').write_text('\n'.join(md)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['results','sample','out']:p.add_argument('--'+key,required=True)
    p.add_argument('--sites',type=int,default=500)
    a=p.parse_args();run(a.results,a.sample,a.out,a.sites)
