"""Construct masked diploid consensus without converting missing sites to reference."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path


def bed(path, lengths):
    result={n:[] for n in lengths}
    if not path or Path(path).name.startswith('NO_'):
        return None
    with open(path,encoding='utf-8') as f:
        for line in f:
            if not line.strip() or line.startswith(('#','track','browser')):continue
            fields=line.split()
            if len(fields)<3 or fields[0] not in lengths:
                raise ValueError('BED contig absent from exact PSMC reference: '+line.strip())
            lo,hi=map(int,fields[1:3]);name=fields[0]
            if not 0<=lo<hi<=lengths[name]:raise ValueError('Out-of-range PSMC BED interval')
            result[name].append((lo,hi))
    for name,rows in result.items():
        merged=[]
        for lo,hi in sorted(rows):
            if merged and lo<=merged[-1][1]:merged[-1]=(merged[-1][0],max(hi,merged[-1][1]))
            else:merged.append((lo,hi))
        result[name]=merged
    return result


def runs(mask):
    import numpy as np
    changes=np.flatnonzero(np.diff(np.concatenate(([0],mask.astype('int8'),[0]))))
    return zip(changes[::2].tolist(),changes[1::2].tolist())


IUPAC={frozenset('AC'):'M',frozenset('AG'):'R',frozenset('AT'):'W',
       frozenset('CG'):'S',frozenset('CT'):'Y',frozenset('GT'):'K'}


def genotype(alleles, gt):
    if gt is None or len(gt)!=2 or any(x is None or x<0 or x>=len(alleles) for x in gt):return None
    bases=[alleles[x].upper() for x in gt]
    if any(len(x)!=1 or x not in 'ACGT' for x in bases):return None
    return bases[0] if bases[0]==bases[1] else IUPAC[frozenset(bases)]


def genotype_confidence(sample):
    gq=sample.get('GQ')
    if gq is not None:return float(gq),'GQ'
    # bcftools does not emit GQ at many monomorphic sites. Retained unseen-allele
    # PLs permit a likelihood-separation guard; this is explicitly not posterior GQ.
    gt=sample.get('GT');pl=sample.get('PL')
    if gt is None or len(gt)!=2 or None in gt or pl is None or any(x is None for x in pl):return None,'unavailable'
    i,j=sorted(gt);called=j*(j+1)//2+i
    if called>=len(pl) or len(pl)<2:return None,'unavailable'
    return float(min(x for k,x in enumerate(pl) if k!=called)-pl[called]),'PL_separation'


def prepare(a):
    import numpy as np
    import pysam
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    digest=hashlib.sha256()
    with open(a.reference,'rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):digest.update(chunk)
    if digest.hexdigest()!=a.reference_sha256:raise ValueError('PSMC reference checksum mismatch')
    fa=pysam.FastaFile(a.reference);lengths=dict(zip(fa.references,fa.lengths))
    inclusion=bed(a.include,lengths);exclusion=bed(a.exclude,lengths);confident=bed(a.confident,lengths)
    if a.mode=='assemblies' and confident is None:raise ValueError('Dipcall confident BED is required')
    vcf=pysam.VariantFile(a.vcf)
    if len(vcf.header.samples)!=1:raise ValueError('PSMC requires one diploid VCF sample')
    for n in vcf.header.contigs:
        if n not in lengths:raise ValueError('VCF/reference contig mismatch: '+n)
        length=vcf.header.contigs[n].length
        if length is not None and length!=lengths[n]:raise ValueError('VCF/reference length mismatch')
    mean_depth=None
    if a.mode=='hifi':
        with open(a.depth) as f:depths=list(csv.DictReader(f,delimiter='\t'))
        selected=[r for r in depths if int(r['endpos'])>=a.min_contig_bp]
        size=sum(int(r['endpos'])-int(r['startpos'])+1 for r in selected)
        mean_depth=sum((int(r['endpos'])-int(r['startpos'])+1)*float(r['meandepth']) for r in selected)/size if size else 0
        minimum=a.min_depth or max(1,math.floor(mean_depth/3))
        maximum=a.max_depth or max(minimum,math.ceil(mean_depth*2))
        if maximum<minimum:raise ValueError('Maximum PSMC depth is below minimum')
    else:minimum=maximum=None
    rows=[];audit=[];confidence_counts={'GQ':0,'PL_separation':0,'unavailable':0};overlapping=0
    with gzip.open(out/'diploid.fa.gz','wt',encoding='ascii') as consensus, (out/'callable.bed').open('w') as maskout:
        for name,length in lengths.items():
            if length<a.min_contig_bp:
                rows.append(dict(contig=name,length_bp=length,scope_bp=0,callable_bp=0,heterozygous_snps=0,status='below_min_contig_bp'))
                continue
            ref=np.frombuffer(fa.fetch(name).upper().encode(),dtype='S1')
            scope=np.zeros(length,dtype=bool)
            for lo,hi in (inclusion[name] if inclusion is not None else [(0,length)]):scope[lo:hi]=True
            if exclusion is not None:
                for lo,hi in exclusion[name]:scope[lo:hi]=False
            scope &= np.isin(ref,[b'A',b'C',b'G',b'T'])
            scope_bp=int(scope.sum())
            seq=np.full(length,b'N',dtype='S1')
            if confident is not None:
                conf=np.zeros(length,dtype=bool)
                for lo,hi in confident[name]:conf[lo:hi]=True
                seq[scope & conf]=ref[scope & conf]
            baseline=int((seq!=b'N').sum())
            bad=[];duplicates=[];snps=0;previous=-1;record_counts={'complex':0,'duplicate':0,'simple':0}
            records=vcf.fetch(name) if name in vcf.header.contigs else []
            for r in records:
                pos=r.start
                if pos<0 or pos+len(r.ref)>length:raise ValueError('VCF record outside reference')
                if ref[pos:pos+len(r.ref)].tobytes().decode()!=r.ref.upper():raise ValueError('VCF REF disagrees with final reference')
                if pos<previous:raise ValueError('Unsorted VCF records')
                if pos==previous:
                    duplicates.append((max(0,pos-a.indel_padding),min(length,r.stop+a.indel_padding)))
                    record_counts['duplicate']+=1
                    overlapping+=1
                    continue
                previous=pos
                sample=r.samples[0];gt=sample.get('GT');alleles=r.alleles or ()
                # Even filtered complex records exclude their reference span and flanks.
                complex_site=any((len(x)!=1 or x.upper() not in 'ACGT') and x not in ['<*>','<NON_REF>'] for x in alleles)
                if complex_site:
                    bad.append((max(0,pos-a.indel_padding),min(length,r.stop+a.indel_padding)))
                    record_counts['complex']+=1
                    continue
                record_counts['simple']+=1
                passed=not (set(r.filter)-{'PASS','.'})
                base=genotype(alleles,gt) if passed else None
                if a.mode=='hifi':
                    dp=sample.get('DP');gq,confidence_source=genotype_confidence(sample)
                    confidence_counts[confidence_source]+=1
                    ad=sample.get('AD')
                    if dp is None or gq is None or not minimum<=dp<=maximum or gq<a.min_gq:base=None
                    if base and base not in 'ACGT':
                        if ad is None or any(x is None for x in ad) or not sum(ad):base=None
                        elif not a.min_allele_balance<=ad[gt[0]]/sum(ad)<=1-a.min_allele_balance:base=None
                # A variant outside Dipcall's confident BED must remain unknown.
                allowed=scope[pos] and (confident is None or conf[pos])
                seq[pos]=(base or 'N').encode() if allowed else b'N'
            before_complex=int((seq!=b'N').sum())
            for lo,hi in bad:seq[lo:hi]=b'N'
            after_complex=int((seq!=b'N').sum())
            for lo,hi in duplicates:seq[lo:hi]=b'N'
            callable_mask=seq!=b'N'
            # Audit the emitted consensus at base resolution, before binning.
            outside_confident_callable_bp=int(np.count_nonzero(callable_mask & ~conf)) if confident is not None else ''
            for lo,hi in runs(callable_mask):maskout.write(f'{name}\t{lo}\t{hi}\n')
            callable_bp=int(callable_mask.sum());snps=int(np.isin(seq,[x.encode() for x in 'MRWSYK']).sum())
            rows.append(dict(contig=name,length_bp=length,scope_bp=int(scope.sum()),callable_bp=callable_bp,
                             heterozygous_snps=snps,status='assessed'))
            audit.append(dict(contig=name,length_bp=length,scope_bp=scope_bp,
                confident_scope_bp=baseline if confident is not None else '',
                outside_confident_bp=scope_bp-baseline if confident is not None else '',
                outside_confident_callable_bp=outside_confident_callable_bp,
                after_simple_genotypes_bp=before_complex,
                simple_genotype_loss_bp=baseline-before_complex if confident is not None else '',
                complex_mask_loss_bp=before_complex-after_complex,
                duplicate_mask_loss_bp=after_complex-callable_bp,
                callable_bp=callable_bp,heterozygous_snps=snps,
                simple_records=record_counts['simple'],complex_records=record_counts['complex'],
                duplicate_records=record_counts['duplicate']))
            consensus.write('>'+name+'\n')
            data=seq.tobytes().decode()
            for i in range(0,length,80):consensus.write(data[i:i+80]+'\n')
    with (out/'contig_qc.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['contig'],delimiter='\t');w.writeheader();w.writerows(rows)
    with (out/'mask_audit.tsv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(audit[0]) if audit else ['contig'],delimiter='\t');w.writeheader();w.writerows(audit)
    total=sum(r['scope_bp'] for r in rows);called=sum(r['callable_bp'] for r in rows);hets=sum(r['heterozygous_snps'] for r in rows)
    status=dict(mode=a.mode,reference_sha256=digest.hexdigest(),coordinate_stage='final',scope_bp=total,
                callable_bp=called,callable_fraction=called/total if total else 0,heterozygous_snps=hets,
                heterozygosity=hets/called if called else None,mean_mapping_depth=mean_depth,
                min_depth=minimum,max_depth=maximum,min_contig_bp=a.min_contig_bp,
                min_gq=a.min_gq,genotype_confidence_counts=confidence_counts,
                duplicate_position_records_masked=overlapping,
                mask_audit_totals={k:sum(r[k] for r in audit) for k in audit[0]
                    if k!='contig' and all(isinstance(r[k],int) for r in audit)} if audit else {},
                genotype_confidence_rule='GQ when provided; otherwise PL runner-up minus called-genotype PL (likelihood separation, not posterior GQ)',
                min_allele_balance=a.min_allele_balance,indel_padding=a.indel_padding,
                include_bed=a.include,exclude_bed=a.exclude,
                scope_note='Long reference scaffolds; sex-linked sequence may be included; no autosomal identity inferred')
    (out/'input_qc.json').write_text(json.dumps(status,indent=2)+'\n')


def finalize(a):
    out=Path(a.out);qc=json.loads((out/'input_qc.json').read_text())
    counts={'K':0,'T':0,'N':0};segments=0;retained=set()
    for line in (out/'input.psmcfa').read_text().splitlines():
        if line.startswith('>'):
            name=line[1:].split()[0]
            if name in retained:raise ValueError('Duplicate PSMC sequence')
            retained.add(name);segments+=1;continue
        if set(line)-set(counts):raise ValueError('Invalid fq2psmcfa output alphabet')
        for c in counts:counts[c]+=line.count(c)
    usable=counts['K']+counts['T']
    with (out/'contig_qc.tsv').open(encoding='utf-8') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    if retained-set(r['contig'] for r in rows):raise ValueError('PSMC sequence absent from input scope')
    retained_callable=sum(int(r['callable_bp']) for r in rows if r['contig'] in retained)
    retained_hets=sum(int(r.get('heterozygous_snps',0)) for r in rows if r['contig'] in retained)
    conversion=[]
    for r in rows:
        keep=r['contig'] in retained
        length=int(r.get('length_bp',0));called=int(r['callable_bp'])
        reason='retained' if keep else ('below_min_contig_bp' if r.get('status')=='below_min_contig_bp'
            else 'below_20_percent_callable' if length and called/length<.2 else 'below_min_good_bases')
        conversion.append(dict(contig=r['contig'],length_bp=length,callable_bp=called,
            callable_fraction=called/length if length else '',retained='yes' if keep else 'no',reason=reason))
    with (out/'conversion_audit.tsv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['contig','length_bp','callable_bp','callable_fraction','retained','reason'],delimiter='\t')
        w.writeheader();w.writerows(conversion)
    reasons=[]
    if qc['callable_bp']<a.min_callable_bp:reasons.append('insufficient_callable_bp')
    if retained_callable<a.min_callable_bp:reasons.append('insufficient_callable_bp_retained_by_fq2psmcfa')
    if qc['callable_fraction']<a.min_callable_fraction:reasons.append('insufficient_callable_fraction')
    if qc['heterozygous_snps']<a.min_hets:reasons.append('insufficient_heterozygous_snps')
    if not counts['K'] or not counts['T']:reasons.append('PSMC_requires_homozygous_and_heterozygous_bins')
    qc.update(status='ready' if not reasons else 'not_fitted',reason='; '.join(reasons),bin_size=a.bin_size,
              bins=counts,usable_bins=usable,segments=segments,psmc_retained_callable_bp=retained_callable,
              retained_contigs=sorted(retained),psmc_retained_heterozygous_snps=retained_hets,
              conversion_lost_callable_bp=qc['callable_bp']-retained_callable,
              fq2psmcfa_min_sequence_callable_fraction=.2,
              missing_bin_fraction=counts['N']/sum(counts.values()) if sum(counts.values()) else 1)
    (out/'input_qc.json').write_text(json.dumps(qc,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare')
    for name in ['reference','reference-sha256','vcf','mode','out']:q.add_argument('--'+name,required=True)
    for name in ['confident','include','exclude','depth']:q.add_argument('--'+name)
    for name,default in [('min-contig-bp',1000000),('min-depth',0),('max-depth',0),('min-gq',20),('indel-padding',5)]:q.add_argument('--'+name,type=int,default=default)
    q.add_argument('--min-allele-balance',type=float,default=.2)
    q=sub.add_parser('finalize');q.add_argument('--out',required=True)
    for name,default in [('bin-size',100),('min-callable-bp',10000000),('min-hets',1000)]:q.add_argument('--'+name,type=int,default=default)
    q.add_argument('--min-callable-fraction',type=float,default=.1)
    a=p.parse_args();prepare(a) if a.command=='prepare' else finalize(a)


if __name__=='__main__':main()
