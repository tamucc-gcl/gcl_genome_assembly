"""Removable comparison of actual PSMC input genotypes and sampled HiFi support."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import subprocess
import sys
import tarfile

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'py_scripts'))
from psmc_input import bed

ALLELES={c:set(c) for c in 'ACGT'}
ALLELES.update(M=set('AC'),R=set('AG'),W=set('AT'),S=set('CG'),Y=set('CT'),K=set('GT'))


def fasta(path):
    name=None;parts=[]
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:yield name,''.join(parts).upper().encode()
                name=line[1:].split()[0];parts=[]
            else:parts.append(line.strip())
    if name is not None:yield name,''.join(parts).upper().encode()


def write_table(path,rows,fields):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)


def support(code,counts,minimum=10):
    if code not in ALLELES:return 'not_called'
    total=sum(counts.values())
    if total<minimum:return 'insufficient_reads'
    alleles=ALLELES[code]
    if sum(counts.get(b,0) for b in alleles)/total<.9:return 'not_supported'
    if len(alleles)==2 and min(counts.get(b,0)/total for b in alleles)<.2:return 'not_supported'
    return 'supported'


def read_support(bam,sites,lengths,genotype_columns=('left','right')):
    """Direct molecule counts shared by paired and single-input diagnostic audits."""
    import pysam
    with pysam.AlignmentFile(str(bam)) as reads:
        if dict(zip(reads.references,reads.lengths))!=lengths:raise ValueError('HiFi BAM sequence dictionary differs')
        for row in sites:
            name=row['scaffold'];pos=row['pos_1based']-1;molecules={}
            for col in reads.pileup(name,pos,pos+1,truncate=True,stepper='samtools',min_mapping_quality=20,min_base_quality=20,max_depth=10000,compute_baq=False):
                if col.reference_pos!=pos:continue
                for item in col.pileups:
                    if item.is_del or item.is_refskip or item.query_position is None:continue
                    read=item.alignment
                    if read.mapping_quality<20 or read.query_qualities is None or read.query_qualities[item.query_position]<20:continue
                    value=read.query_sequence[item.query_position].upper()
                    if value in 'ACGT':molecules.setdefault(read.query_name,set()).add(value)
            bases=Counter(next(iter(v)) for v in molecules.values() if len(v)==1)
            result=dict(row,molecules=sum(bases.values()),**{b:bases[b] for b in 'ACGT'})
            result.update({key+'_support':support(row[key],bases) for key in genotype_columns})
            yield result


def hifi_bam(run,sample,reference_hash):
    run=Path(run);plan=json.loads((run/'results/post_assembly/psmc/psmc_plan.json').read_text())
    job=next(j for j in plan['jobs'] if j['sample']==sample and j['mode']=='hifi')
    if job['reference']['sha256']!=reference_hash:raise ValueError('HiFi comparison uses a different reference')
    with (run/'trace.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    hits=[r for r in rows if r['name'].split(' (')[0].endswith(':PSMC_MAP_HIFI') and r['name'].endswith('('+sample+')')
          and r['status'] in ['COMPLETED','CACHED'] and r['exit']=='0']
    if len(hits)!=1:raise ValueError('Expected one successful same-individual HiFi mapping in supplied trace')
    root=Path(hits[0]['workdir']);bam=root/'hifi.bam'
    if not bam.is_file() or not Path(str(bam)+'.bai').is_file():raise ValueError('Recorded HiFi BAM/index unavailable')
    digest=hashlib.sha256()
    with (root/'reference.fa').open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):digest.update(chunk)
    if digest.hexdigest()!=reference_hash:raise ValueError('Mapped reference content no longer matches')
    return bam


def compare_block(a,b,scope):
    import numpy as np
    ca=np.isin(a,[c.encode() for c in ALLELES]);cb=np.isin(b,[c.encode() for c in ALLELES])
    ha=np.isin(a,[c.encode() for c in 'MRWSYK']);hb=np.isin(b,[c.encode() for c in 'MRWSYK'])
    masks=dict(shared=scope & ca & cb,left_only=scope & ca & ~cb,
        right_only=scope & ~ca & cb,neither=scope & ~ca & ~cb)
    rows=[]
    for label,mask in masks.items():
        rows.append(dict(region=label,bp=int(mask.sum()),left_hets=int((mask & ha).sum()),right_hets=int((mask & hb).sum()),
            discordant=int((mask & (a!=b) & ca & cb).sum()),het_vs_hom=int((mask & ca & cb & (ha!=hb)).sum()),
            either_het=int((mask & (ha|hb)).sum()),discordant_with_het=int((mask & ca & cb & (a!=b) & (ha|hb)).sum())))
    candidates=dict(shared_discordant=masks['shared'] & (a!=b),shared_concordant_het=masks['shared'] & (a==b) & ha,
        left_only_het=masks['left_only'] & ha,right_only_het=masks['right_only'] & hb)
    return rows,candidates


def main():
    import numpy as np
    import pysam
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--comparison',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--hifi-run',type=Path);p.add_argument('--sample-sites',type=int,default=500)
    p.add_argument('--left-run',default='asm10');p.add_argument('--right-run',default='asm20')
    a=p.parse_args()
    if a.sample_sites<1:raise ValueError('Positive site sample size required')
    if a.left_run==a.right_run:raise ValueError('Select two distinct comparison runs')
    root=a.comparison.resolve();out=a.out.resolve()
    if out.exists() and any(out.iterdir()):raise ValueError('Use a new reliability output directory')
    comparison=json.loads((root/'comparison.json').read_text());sample=comparison['sample']
    source=json.loads((Path(comparison['source'])/'manifest.json').read_text());job=source['controls'][0]
    dirs=[];qcs=[];provenances=[]
    for preset in [a.left_run,a.right_run]:
        run=Path(comparison.get('run_directories',{}).get(preset,str(root/preset)))
        if (run/'exit_status.txt').read_text().strip()!='0':raise ValueError('Full comparison did not succeed: '+preset)
        plan=json.loads((run/'results/post_assembly/psmc/psmc_plan.json').read_text())
        chosen=next(j for j in plan['jobs'] if j['sample']==sample)
        base=run/'results/post_assembly/psmc'/str(chosen['taxid'])/sample;dirs.append(base)
        qcs.append(json.loads((base/'input/input_qc.json').read_text()))
        provenances.append(json.loads((base/'input/provenance.json').read_text()))
    if qcs[0]['reference_sha256']!=qcs[1]['reference_sha256'] or provenances[0]['filters']!=provenances[1]['filters']:
        raise ValueError('References or input filters differ')
    if provenances[0]['filters']!=job['filters'] or job['reference']['sha256']!=qcs[0]['reference_sha256']:
        raise ValueError('Full comparisons differ from audited mask scope/settings')
    bam=hifi_bam(a.hifi_run,sample,qcs[0]['reference_sha256']) if a.hifi_run else None
    out.mkdir(parents=True,exist_ok=True);totals={k:Counter() for k in ['shared','left_only','right_only','neither']}
    windows=[];reservoir={};rng=np.random.default_rng(1729)
    reference=job['reference']['fasta']
    digest=hashlib.sha256()
    with open(reference,'rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):digest.update(chunk)
    if digest.hexdigest()!=qcs[0]['reference_sha256']:raise ValueError('Audited final reference changed')
    with pysam.FastaFile(reference) as fa:
        lengths=dict(zip(fa.references,fa.lengths));inc=bed(job.get('include_bed'),lengths);exc=bed(job.get('exclude_bed'),lengths)
        for left,right in zip_longest(fasta(dirs[0]/'input/diploid.fa.gz'),fasta(dirs[1]/'input/diploid.fa.gz')):
            if left is None or right is None or left[0]!=right[0] or len(left[1])!=len(right[1]):raise ValueError('Consensus sequence dictionaries differ')
            name,seqA=left;seqB=right[1]
            if len(seqA)!=lengths[name]:raise ValueError('Consensus/reference coordinate mismatch')
            for lo in range(0,len(seqA),1000000):
                hi=min(len(seqA),lo+1000000);x=np.frombuffer(seqA[lo:hi],dtype='S1');y=np.frombuffer(seqB[lo:hi],dtype='S1')
                ref=np.frombuffer(fa.fetch(name,lo,hi).upper().encode(),dtype='S1');scope=np.isin(ref,[b'A',b'C',b'G',b'T'])
                if inc is not None:
                    included=np.zeros(hi-lo,dtype=bool)
                    for start,end in inc[name]:
                        if end>lo and start<hi:included[max(0,start-lo):min(hi,end)-lo]=True
                    scope &= included
                if exc is not None:
                    for start,end in exc[name]:
                        if end>lo and start<hi:scope[max(0,start-lo):min(hi,end)-lo]=False
                rows,candidates=compare_block(x,y,scope)
                for r in rows:
                    totals[r['region']].update({k:v for k,v in r.items() if k!='region'})
                    windows.append(dict(scaffold=name,start=lo,end=hi,**r))
                for label,mask in candidates.items():
                    positions=np.flatnonzero(mask);keys=rng.random(len(positions))
                    if len(keys)>a.sample_sites:
                        selected=np.argpartition(keys,a.sample_sites-1)[:a.sample_sites];keys=keys[selected];positions=positions[selected]
                    values=reservoir.get(label,[])+[(float(key),name,int(lo+pos),ref[pos].decode(),x[pos].decode(),y[pos].decode()) for key,pos in zip(keys,positions)]
                    reservoir[label]=sorted(values)[:a.sample_sites]
    shared=totals['shared']['bp']
    if shared+totals['left_only']['bp']!=qcs[0]['callable_bp'] or shared+totals['right_only']['bp']!=qcs[1]['callable_bp']:
        raise ValueError('Consensus comparison does not reproduce final callable counts')
    if sum(r['bp'] for r in totals.values())!=qcs[0]['scope_bp']:raise ValueError('Scope count mismatch')
    counts=[dict(region=k,**v) for k,v in totals.items()]
    fields=['region','bp','left_hets','right_hets','discordant','het_vs_hom','either_het','discordant_with_het']
    write_table(out/'regions.tsv',counts,fields);write_table(out/'windows.tsv',windows,['scaffold','start','end']+fields)
    sites=[]
    if bam:
        selected=[dict(group=group,scaffold=name,pos_1based=pos+1,reference=ref,left=x,right=y)
            for group,items in reservoir.items() for _,name,pos,ref,x,y in sorted(items,key=lambda r:(r[1],r[2]))]
        sites=list(read_support(bam,selected,lengths))
        write_table(out/'sampled_read_support.tsv',sites,['group','scaffold','pos_1based','reference','left','right','molecules','A','C','G','T','left_support','right_support'])
    else:
        for group,items in reservoir.items():
            sites.extend(dict(group=group,scaffold=n,pos_1based=pos+1,reference=r,left=x,right=y) for _,n,pos,r,x,y in items)
        write_table(out/'sampled_sites.tsv',sites,['group','scaffold','pos_1based','reference','left','right'])
    summary=dict(schema_version=2,sample=sample,labels=dict(left=a.left_run,right=a.right_run),regions=counts,scope_bp=qcs[0]['scope_bp'],shared_genotype_agreement=1-totals['shared']['discordant']/shared,
        shared_het_genotype_agreement=1-totals['shared']['discordant_with_het']/totals['shared']['either_het'] if totals['shared']['either_het'] else None,
        reference_sha256=qcs[0]['reference_sha256'],hifi_bam=str(bam) if bam else None,seed=1729,sites_per_group=a.sample_sites,
        scope='Actual masked diploid consensus before PSMC binning',retained_scaffolds=[q['retained_contigs'] for q in qcs])
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    md=['# PSMC genotype reliability comparison','',f'Individual: **{sample}**. Left: **{a.left_run}**. Right: **{a.right_run}**.','',
        f'Shared callable genotype agreement: **{100*summary["shared_genotype_agreement"]:.4f}%**.','',
        f'| Region | Mb | {a.left_run} heterozygosity | {a.right_run} heterozygosity | Discordant calls |','| --- | ---: | ---: | ---: | ---: |']
    region_labels=dict(shared='Callable in both',left_only='Only '+a.left_run,right_only='Only '+a.right_run,neither='Callable in neither')
    for r in counts:
        rate=lambda key:f"{100*r[key]/r['bp']:.3f}%" if r['bp'] else 'NA'
        md.append(f"| {region_labels[r['region']]} | {r['bp']/1e6:.2f} | {rate('left_hets') if r['region'] in ['shared','left_only'] else 'NA'} | {rate('right_hets') if r['region'] in ['shared','right_only'] else 'NA'} | {r['discordant']} |")
    md+=['','![Shared-scope genotype comparison](reliability.png)','',
        ('Agreement at shared callable positions heterozygous in either preset: **'+f'{100*summary["shared_het_genotype_agreement"]:.3f}%**.' if summary['shared_het_genotype_agreement'] is not None else 'No shared callable heterozygous positions.'),
        'This heterozygote-focused measure prevents abundant agreeing homozygous-reference bases from hiding differences relevant to PSMC.','',
        'The comparison uses the actual unphased IUPAC genotypes entering PSMC, including homozygous reference calls. Counts independently reproduce both final callable masks. Aggregate differences can arise from both scope selection and changed genotypes.','',
        '[Region counts](regions.tsv) · [1 Mb windows](windows.tsv) · [Provenance](summary.json)','',
        'Read support, when supplied, is a diagnostic on deterministic random samples, not a truth set or an automatic preset choice. It uses primary/supplementary placements passing MAPQ 20 and base quality 20, default pileup exclusion flags, and distinct read names; conflicting bases from the same molecule are excluded. BAQ is disabled for this direct base-count diagnostic, unlike the production caller. At least ten observed molecules, ≥90% bases matching the proposed alleles and ≥20% support for each heterozygous allele define the displayed support label.','',
        'The BAM uses a single haplotype reference, so mapping bias and missing alternate placements remain limitations. The same HiFi data also contributed to assembly; these are direct read checks, not independent sequencing validation.']
    if bam:
        md+=['',f'| Sample group | Sites | {a.left_run} supported | {a.right_run} supported | Both supported | Neither supported with ≥10 molecules | Insufficient molecules |',
             '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
        for group in reservoir:
            rs=[s for s in sites if s['group']==group]
            md.append(f"| {group} | {len(rs)} | {sum(s['left_support']=='supported' for s in rs)} | {sum(s['right_support']=='supported' for s in rs)} | {sum(s['left_support']==s['right_support']=='supported' for s in rs)} | {sum(s['molecules']>=10 and 'supported' not in [s['left_support'],s['right_support']] for s in rs)} | {sum(s['molecules']<10 for s in rs)} |")
        md+=['','[Sampled read evidence](sampled_read_support.tsv)']
    else:md+=['','No HiFi BAM was supplied; genotype concordance alone cannot establish which discordant call is correct. [Sampled sites](sampled_sites.tsv)']
    (out/'README.md').write_text('\n'.join(md)+'\n')
    print('Report:',out/'README.md')


if __name__=='__main__':main()
