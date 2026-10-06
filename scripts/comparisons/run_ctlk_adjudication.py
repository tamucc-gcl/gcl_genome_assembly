#!/usr/bin/env python3
"""Read-only, parallel junction evidence lanes. No automatic sequence cuts."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import csv
import re
import shutil
import subprocess

from trace_chimera_origins import Fasta, run
from chimera_origin import table, sha, agp_rows, revcomp
from assess_junction_batch import write_rows, pairs_dictionary, verify_agp
from junction_assessment import Lift, library_for
from junction_focus import read_blocks, intersection_length

CORES = {
    1: [('J01', 'h1tg000004l', 109308, 433057),
        ('J02', 'h1tg000133l', 4022347, 4192059),
        ('J05', 'h1tg000153l', 719979, 857651)],
    2: [('J06', 'h2tg000028l', 660967, 839165),
        ('J07', 'h2tg000298l', 638671, 777533)]}


def paf(path):
    with Path(path).open() as src:
        for line in src:
            f = line.rstrip().split('\t')
            if len(f) < 12:
                raise ValueError('Malformed PAF')
            yield dict(query=f[0], query_length=int(f[1]), qs=int(f[2]), qe=int(f[3]),
                       strand=f[4], target=f[5], target_length=int(f[6]),
                       ts=int(f[7]), te=int(f[8]), matches=int(f[9]), block=int(f[10]),
                       mapq=int(f[11]), tags=';'.join(f[12:]))


def classify_anchor(hits, length):
    """Conservative placement screen; MAPQ alone is insufficient for repeats."""
    good = [h for h in hits if h['qe']-h['qs'] >= .9*length and h['matches']/h['block'] >= .9]
    good.sort(key=lambda h: h['matches'], reverse=True)
    if not good:
        return 'no_full_placement', None
    best = good[0]
    competitors = [h for h in good[1:] if (h['target'], h['strand']) != (best['target'], best['strand'])
                   or abs(h['ts']-best['ts']) > length]
    if best['mapq'] < 30 or any(h['matches'] >= .95*best['matches'] for h in competitors):
        return 'ambiguous', best
    return 'unique_screen', best


def map_fasta(a, target, query, prefix, reads=False):
    cmd = ['minimap2', '-t', str(a.threads), '-I', '8G', '-x', 'map-hifi' if reads else 'asm5',
           '--secondary=yes', '-N', '50', '-p', '0.5', '-c', '--eqx', str(target), str(query)]
    run(cmd, stdout=str(prefix)+'.paf', stderr=str(prefix)+'.log')
    return list(paf(str(prefix)+'.paf'))


def peers(a):
    """Use current sisters and six better haplotypes, never historical assemblies."""
    paths = []
    for sample in ('Sde-CBau_104', 'Sde-CLim_110', 'Sde-CMat_203', 'Sde-CTlk_104'):
        for h in (1, 2):
            p = a.assessment/'assembly/scaffold/yahs_round2'/f'{sample}_hap{h}_round2_scaffolds.fa'
            if not p.is_file():
                raise ValueError(f'Missing current comparison FASTA: {p}')
            paths.append(p)
    return paths


def setup(a):
    status = json.loads((a.packet/'status.json').read_text())
    if status['status'] != 'SUCCESS':
        raise ValueError('Incomplete source evidence packet')
    job = next(j for j in status['assemblies'] if j['assembly'].endswith(f'hap{a.hap}'))
    folder = a.packet/job['folder']
    origin = json.loads((folder/'origin_provenance.json').read_text())
    final = next(s for s in origin['stages'] if s['stage'] == 'assessment')
    if sha(Path(final['fasta'])) != job['assessment_sha256']:
        raise ValueError('Assessment FASTA changed')
    fa = Fasta(Path(final['fasta']), a.out/'assessment.fa', 'samtools')
    regions = []
    # Exact existing candidate intervals and their verified gaps/controls.
    for name, r in origin['intervals'].items():
        regions.append(dict(id=name, scaffold=r['scaffold'], start=int(r['lo']), end=int(r['hi']), role='transition'))
    for r in table(folder/'selected_gaps.tsv'):
        regions.append(dict(id=r['id'], scaffold=r['scaffold'], start=int(r['start']),
                            end=int(r['end']), role=r['role']))
    for r in regions:
        if not 0 <= r['start'] < r['end'] <= fa.lengths[r['scaffold']]:
            raise ValueError('Region outside exact assessment')
    write_rows(a.out/'regions.tsv', regions)
    shutil.copy2(folder/'origin_provenance.json', a.out/'origin_provenance.json')
    return job, folder, final, fa, regions


def anchors(a, fa, regions):
    metadata = []
    with (a.out/'anchors.fa').open('w') as dest:
        def emit(chrom, lo, hi, region, side):
            if not 0 <= lo < hi <= fa.lengths[chrom]:
                return
            seq = fa.fetch(chrom, lo, hi)
            q = f'a{len(metadata):06d}'
            metadata.append(dict(query=q, region=region, side=side, scaffold=chrom,
                                 start=lo, end=hi, length=hi-lo, ambiguous_bases=sum(b not in 'ACGT' for b in seq.upper())))
            dest.write(f'>{q}\n{seq}\n')
        for r in regions:
            for size in (10000, 25000):
                for offset in (0, 25000, 100000, 250000, 500000, 1000000):
                    emit(r['scaffold'], r['start']-offset-size, r['start']-offset, r['id'], 'left')
                    emit(r['scaffold'], r['end']+offset, r['end']+offset+size, r['id'], 'right')
        # Broad chromosome-piece correspondence, including suffix pieces.
        for chrom, length in fa.lengths.items():
            if length >= 5000000:
                for lo in range(0, length-25000+1, 1000000):
                    emit(chrom, lo, lo+25000, 'chromosome_piece', 'tile')
    write_rows(a.out/'anchor_coordinates.tsv', metadata)
    byquery = {r['query']: r for r in metadata}
    summary = []
    for i, peer in enumerate(peers(a)):
        print(f'Anchors: {peer.name}', flush=True)
        hits = map_fasta(a, peer, a.out/'anchors.fa', a.out/f'peer_{i}')
        grouped = defaultdict(list)
        for h in hits:
            grouped[h['query']].append(h)
        # Name maps are copied for exact chromosome-label interpretation.
        name = peer.name.split('_round2')[0]
        nm = a.assessment/'assembly/harmonization'/f'{name}.harmonized_name_map.tsv'
        if not nm.is_file():
            raise ValueError(f'Missing chromosome name map: {nm}')
        shutil.copy2(nm, a.out/f'peer_{i}.name_map.tsv')
        labels = defaultdict(set)
        for r in table(nm):
            labels[r['old_name']].add(r['new_name'])
        for q, meta in byquery.items():
            state, h = classify_anchor(grouped[q], meta['length'])
            if meta['ambiguous_bases']:
                state = 'ambiguous_sequence'
            summary.append(dict(**meta, peer=peer.name, individual=name.rsplit('_hap', 1)[0],
                                status=state, target=h['target'] if h else '.',
                                chromosome_label=';'.join(sorted(labels[h['target']])) if h else '.',
                                target_start=h['ts'] if h else '.', target_end=h['te'] if h else '.',
                                strand=h['strand'] if h else '.', mapq=h['mapq'] if h else '.',
                                emitted_full_placements=len(grouped[q])))
    write_rows(a.out/'anchor_correspondence.tsv', summary)
    return {'anchors': len(metadata), 'interpretation': 'Conservative screen, not formal uniqueness; emitted alternatives capped at 50. Inspect PAF competitors and consistent independent tiles.'}


def graph_neighborhood(path, seeds):
    """One-hop undirected discovery; exported native links keep their orientations."""
    names = set(seeds)
    with Path(path).open() as src:
        for line in src:
            f = line.rstrip().split('\t')
            if f[0] == 'L' and (f[1] in seeds or f[3] in seeds):
                names.update((f[1], f[3]))
    return names


def graphs(a):
    raw = a.assessment/'assembly/contig/hifiasm'/f'Sde-CTlk_104.hap{a.hap}.p_ctg.fasta'
    fa = Fasta(raw, a.out/'raw.fa', 'samtools')
    with (a.out/'cores.fa').open('w') as dest:
        for j, chrom, lo, hi in CORES[a.hap]:
            dest.write(f'>{j}\n{fa.fetch(chrom, max(0,lo-100000), min(fa.lengths[chrom],hi+100000))}\n')
    summary = []
    for treatment in ('both', 'hifi_only'):
        source = a.fusion/treatment
        if json.loads((source/'status.json').read_text())['status'] != 'SUCCESS':
            raise ValueError('Graph source treatment incomplete')
        paths = sorted(source.glob('*.gfa'))
        if not paths:
            raise ValueError(f'Native graphs missing: {source}')
        for index, path in enumerate(paths):
            print(f'Graph scan: {path.name}', flush=True)
            dest = a.out/f'{treatment}_{index}'; dest.mkdir()
            segmentfa = dest/'segments.fa'
            missing = 0
            with path.open() as src, segmentfa.open('w') as out:
                for line in src:
                    f = line.rstrip().split('\t')
                    if f[0] == 'S':
                        if f[2] == '*':
                            missing += 1
                        else:
                            out.write(f'>{f[1]}\n{f[2]}\n')
            if segmentfa.stat().st_size == 0:
                summary.append(dict(treatment=treatment, graph=str(path), status='no_segment_sequences', missing_sequences=missing))
                continue
            hits = map_fasta(a, segmentfa, a.out/'cores.fa', dest/'core_segments')
            seeds = {h['target'] for h in hits if h['qe']-h['qs'] >= 1000}
            names = graph_neighborhood(path, seeds)
            # Keep small graph topology/tags/read-placement records, omit huge sequences.
            with path.open() as src, (dest/'local.gfa').open('w') as out:
                for line in src:
                    f = line.rstrip().split('\t'); kind=f[0]
                    if kind == 'S' and f[1] in names:
                        tags=f[3:]
                        if f[2]!='*' and not any(t.startswith('LN:i:') for t in tags):
                            tags=tags+[f'LN:i:{len(f[2])}']
                        out.write('\t'.join([f[0],f[1],'*']+tags)+ '\n')
                    elif kind == 'L' and f[1] in names and f[3] in names:
                        out.write(line)
                    elif kind == 'A' and f[1] in names:
                        out.write(line)
                    elif kind in ('P','W'):
                        # Retain paths only if they explicitly mention a selected segment.
                        tokens=(x.rstrip('+-') for x in re.findall(r'[^,><\s\t]+', '\t'.join(f[2:])))
                        if any(x in names for x in tokens):
                            out.write(line)
            write_rows(dest/'seed_placements.tsv', hits)
            summary.append(dict(treatment=treatment, graph=str(path), status='mapped', seeds=len(seeds),
                                neighborhood_segments=len(names), missing_sequences=missing))
            segmentfa.unlink()  # Only this newly created scratch FASTA.
    write_rows(a.out/'graph_inventory.tsv', summary)
    return {'graphs': len(summary), 'interpretation': 'Sequence placement plus native one-hop topology/A records; not an inferred assembly-edge history or uniquely spanning molecule count.'}


def retrieve_sequences(sam):
    """Keep the longest available representation; flag hard-clipped ascertainment."""
    seqs = {}; clipped=set()
    with Path(sam).open() as src:
        for line in src:
            if line.startswith('@'):
                continue
            f=line.rstrip().split('\t')
            if len(f)<11:
                raise ValueError('Malformed SAM')
            if f[9]=='*':
                continue
            if 'H' in f[5]:
                clipped.add(f[0])
            seq=revcomp(f[9]) if int(f[1])&16 else f[9]
            if len(seq)>len(seqs.get(f[0],'')):
                seqs[f[0]]=seq
    return seqs, clipped


def distance_area(first, second, lower, upper):
    """Continuous coordinate-pair exposure for lower <= y-x < upper."""
    a,b=first;c,d=second
    def cdf(t):
        positive=lambda v:max(0,v)**2
        return .5*(positive(t-c+b)-positive(t-c+a)-positive(t-d+b)+positive(t-d+a))
    return max(0,cdf(upper)-cdf(lower))


def hifi(a, job, fa, regions):
    bam=Path(job['hifi']['bam'])
    run(['samtools','view','-H',bam],stdout=a.out/'bam.header.sam')
    dictionary={}
    for line in (a.out/'bam.header.sam').read_text().splitlines():
        if line.startswith('@SQ\t'):
            tags=dict(x.split(':',1) for x in line.split('\t')[1:])
            dictionary[tags['SN']]=int(tags['LN'])
    if dictionary!=fa.lengths:
        raise ValueError('HiFi BAM and assessment reference dictionaries differ')
    # Union BED avoids rereading overlapping loci; view -M avoids repeated records.
    grouped=defaultdict(list)
    for r in regions:
        grouped[r['scaffold']].append((max(0,r['start']-250000),min(fa.lengths[r['scaffold']],r['end']+250000)))
    from junction_focus import union
    with (a.out/'retrieval.bed').open('w') as dest:
        for chrom, intervals in grouped.items():
            for lo,hi in union(intervals):dest.write(f'{chrom}\t{lo}\t{hi}\n')
    run(['samtools','view','-M','-L',a.out/'retrieval.bed',bam],stdout=a.out/'retrieved.sam')
    seqs, clipped=retrieve_sequences(a.out/'retrieved.sam')
    if not seqs:
        raise ValueError('No retrievable HiFi sequences')
    with (a.out/'reads.fa').open('w') as dest:
        for q,seq in seqs.items():dest.write(f'>{q}\n{seq}\n')
    # Candidate windows seed corresponding competing regions; do not create artificial joins.
    windows=[]
    with (a.out/'queries.fa').open('w') as dest:
        for i,r in enumerate(regions):
            lo=max(0,r['start']-100000); hi=min(fa.lengths[r['scaffold']],r['end']+100000)
            q=f'w{i}'; windows.append(dict(query=q,**r,lo=lo,hi=hi))
            dest.write(f'>{q}\n{fa.fetch(r['scaffold'],lo,hi)}\n')
    references=[]; seen=set()
    with (a.out/'competitive.fa').open('w') as dest:
        for i,peer in enumerate(peers(a)):
            pf=Fasta(peer,a.out/f'peer_{i}.fa','samtools')
            placements=map_fasta(a,peer,a.out/'queries.fa',a.out/f'peer_{i}_regions')
            candidates=[]
            if peer.resolve()==Path(fa.path).resolve():
                candidates.extend((r['scaffold'],r['lo'],r['hi']) for r in windows)
            candidates.extend((h['target'],max(0,h['ts']-50000),min(h['target_length'],h['te']+50000))
                              for h in placements if h['qe']-h['qs']>=1000)
            # Merge overlapping windows to avoid artificial duplicates suppressing MAPQ.
            bychrom=defaultdict(list)
            for chrom,lo,hi in candidates:bychrom[chrom].append((lo,hi))
            for chrom,ivals in bychrom.items():
                for lo,hi in union(ivals):
                    key=(str(peer.resolve()),chrom,lo,hi)
                    if key in seen:continue
                    seen.add(key); name=f'c{len(references):05d}'
                    references.append(dict(name=name,assembly=peer.name,scaffold=chrom,start=lo,end=hi))
                    dest.write(f'>{name}\n{pf.fetch(chrom,lo,hi)}\n')
    write_rows(a.out/'competitive_reference.tsv',references)
    placements=map_fasta(a,a.out/'competitive.fa',a.out/'reads.fa',a.out/'read_placements',reads=True)
    write_rows(a.out/'read_placements.tsv',placements)
    # Export aligned read spans at the ORIGINAL exact candidate/control flanks.
    ref_byname={r['name']:r for r in references}; evidence=[]
    for hit in placements:
        ref=ref_byname[hit['target']]
        if ref['assembly']!=Path(fa.path).resolve().name:
            continue
        # PAF target intervals are local; cg blocks exclude deletions when counting anchors.
        cg=next((x[5:] for x in hit['tags'].split(';') if x.startswith('cg:Z:')),None)
        if not cg:raise ValueError('PAF missing CIGAR')
        _,blocks=read_blocks(ref['start']+hit['ts'],cg)
        for r in regions:
            if r['scaffold']!=ref['scaffold']:continue
            left=intersection_length([(x,y) for x,y,_,_ in blocks],[(max(0,r['start']-2000),r['start'])])
            right=intersection_length([(x,y) for x,y,_,_ in blocks],[(r['end'],r['end']+2000)])
            if left or right:
                evidence.append(dict(region=r['id'],role=r['role'],read=hit['query'],mapq=hit['mapq'],
                    left_aligned_bp=left,right_aligned_bp=right,spans_both_1kb=left>=1000 and right>=1000,
                    representation_may_be_clipped=hit['query'] in clipped,target=hit['target'],
                    start=ref['start']+hit['ts'],end=ref['start']+hit['te']))
    write_rows(a.out/'original_flank_read_evidence.tsv',evidence)
    return {'selected_reads':len(seqs),'clipped_names':len(clipped),'placements':len(placements),
            'interpretation':'BAM-ascertained local reads, emitted alternatives capped at 50. Unmapped reads absent from regional retrieval. Longest available sequence can remain hard-clipped. Peer regions are a competitive screen, not phase-equivalent paths. Broad-core zero bridges cannot justify cutting; inspect intact controls and narrower graph boundaries.'}


def hic(a, job, folder, final, fa, regions):
    meta=job['hic']
    if sha(Path(final['agp']))!=meta['agp_sha256'] or sha(Path(meta['source_fasta']))!=meta['source_fasta_sha256']:
        raise ValueError('Hi-C mapping provenance changed')
    source_lengths={f[0]:int(f[1]) for f in (l.split('\t') for l in (folder/'hic_input.fa.fai').read_text().splitlines())}
    if pairs_dictionary(Path(meta['pairs']))!=source_lengths:
        raise ValueError('Hi-C pairs reference dictionary differs')
    if sha(folder/'hic_readsets.tsv')!=meta['readsets_sha256']:
        raise ValueError('Hi-C library manifest changed')
    libs=table(folder/'hic_readsets.tsv'); rows=agp_rows(final['agp']);verify_agp(rows,fa)
    lift=Lift(rows); size=250000
    focal={c for c,n in fa.lengths.items() if n>=5000000}
    counts=Counter(); totals=Counter(); margins=Counter(); scanned=0
    decay=Counter()
    # Sparse region index avoids five Contacts allocations for every genome-wide pair.
    near=defaultdict(list)
    for r in regions:
        lo=max(0,r['start']-1000000);hi=min(fa.lengths[r['scaffold']],r['end']+1000000)
        for b in range(lo//size,(hi-1)//size+1):near[r['scaffold'],b].append((r,lo,hi))
    from junction_assessment import Contacts
    scales={s:Contacts(regions,fa.lengths,s) for s in (25000,100000,250000,500000,1000000)}
    with gzip.open(meta['pairs'],'rt') as src:
        for line in src:
            if line.startswith('#') or not line.strip():continue
            scanned+=1
            if scanned%10000000==0:print(f'Hi-C scanned {scanned:,}',flush=True)
            f=line.rstrip().split('\t')
            if len(f)<8:raise ValueError('Malformed pairs')
            if f[7]!='UU':continue
            x=lift.locate(f[1],int(f[2])-1);y=lift.locate(f[3],int(f[4])-1)
            if x is None or y is None:continue
            lib=library_for(f[0],libs)
            if lib=='UNASSIGNED':raise ValueError('Unassigned Hi-C library')
            totals[lib]+=1
            affected={r['id']:r for z in (x,y) for r,lo,hi in near[z[0],z[1]//size] if lo<=z[1]<hi}
            if affected:
                for c in scales.values():c.add(lib,x,y)
                if x[0]==y[0]:
                    p1,p2=sorted((x[1],y[1]));distance=p2-p1
                    upper=next((b for b in (25000,50000,100000,250000,500000,1000000,2000000) if distance<b),None)
                    if upper:
                        for r in affected.values():
                            if x[0]!=r['scaffold']:continue
                            left=max(0,r['start']-1000000);right=min(fa.lengths[x[0]],r['end']+1000000)
                            if left<=p1<r['start'] and r['end']<=p2<right:kind='cross'
                            elif left<=p1<=p2<r['start']:kind='left_cis'
                            elif r['end']<=p1<=p2<right:kind='right_cis'
                            else:continue
                            decay[lib,r['id'],r['role'],kind,upper]+=1
            if x[0] not in focal and y[0] not in focal:continue
            b1=(x[0],x[1]//size);b2=(y[0],y[1]//size)
            b1,b2=sorted((b1,b2));counts[lib,b1,b2]+=1
            margins[lib,b1]+=1;margins[lib,b2]+=1
    if scanned!=meta['audits']['input_pairs']:raise ValueError('Input pairs count differs from retained packet')
    # Stream the large matrix instead of materializing millions of row dictionaries.
    fields=['library','scaffold1','bin1','scaffold2','bin2','bin_bp','count','marginal1','marginal2','pairs_per_million']
    with gzip.open(a.out/'contact_bins.tsv.gz','wt') as dest:
        writer=csv.DictWriter(dest,fieldnames=fields,delimiter='\t');writer.writeheader()
        for (l,x,y),n in counts.items():
            writer.writerow(dict(library=l,scaffold1=x[0],bin1=x[1],scaffold2=y[0],bin2=y[1],
                bin_bp=size,count=n,marginal1=margins[l,x],marginal2=margins[l,y],pairs_per_million=n*1e6/totals[l]))
    profiles=[]
    bounds=(0,25000,50000,100000,250000,500000,1000000,2000000)
    for lib in sorted(totals):
        for r in regions:
            left=(max(0,r['start']-1000000),r['start'])
            right=(r['end'],min(fa.lengths[r['scaffold']],r['end']+1000000))
            for kind,first,second in (('cross',left,right),('left_cis',left,left),('right_cis',right,right)):
                for lower,upper in zip(bounds,bounds[1:]):
                    n=decay[lib,r['id'],r['role'],kind,upper]
                    area=distance_area(first,second,lower,upper)
                    profiles.append(dict(library=lib,id=r['id'],role=r['role'],kind=kind,
                        distance_bin_lower=lower,distance_bin_upper=upper,count=n,
                        coordinate_pair_area_bp2=area,pairs_per_million=n*1e6/totals[lib],
                        per_million_pairs_per_Mb2=n*1e18/(totals[lib]*area) if area else '.',
                        normalization='distance/window geometry and library size only; not mappability'))
    write_rows(a.out/'distance_profiles.tsv',profiles)
    out=[];partners=[]
    for s,c in scales.items():
        grouped=defaultdict(list)
        for (lib,rid,side,chrom,b),n in c.partners.items():grouped[lib,rid,side].append((n,chrom,b))
        for lib in sorted(totals):
            for r in regions:
                out.append(dict(id=r['id'],role=r['role'],flank_bp=s,library=lib,
                    crossing_pairs=c.cross[lib,r['id']],projected_UU_pairs=totals[lib],
                    left_contacts=c.margin[lib,r['id'],'left'],right_contacts=c.margin[lib,r['id'],'right']))
                for side in ('left','right'):
                    for n,chrom,b in sorted(grouped[lib,r['id'],side],reverse=True)[:20]:
                        partners.append(dict(id=r['id'],side=side,flank_bp=s,library=lib,partner=chrom,start=b*s,count=n))
    write_rows(a.out/'contact_scales.tsv',out);write_rows(a.out/'competing_partners.tsv',partners)
    # Raw bin matrix retains separation for distance/marginal adjustment after unique-anchor annotation.
    return {'input_pairs':scanned,'projected_by_library':dict(totals),'bin_bp':size,
            'interpretation':'Raw library-separated matrix, flanks and margins; no uniqueness filtering beyond upstream UU/MAPQ filtering. Distance/mappability-calibrated adjudication requires anchor results. Reuses scaffolding reads, not independent validation.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--lane',choices=('anchors','graph','hifi','hic'),required=True)
    p.add_argument('--hap',type=int,choices=(1,2),required=True)
    p.add_argument('--packet',type=Path,required=True)
    p.add_argument('--assessment',type=Path,required=True)
    p.add_argument('--fusion',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--threads',type=int,default=8)
    a=p.parse_args()
    if a.threads<=0:p.error('Positive threads required')
    a.out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(Path(__file__),a.out/'runner.py')
    manifest={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}
    manifest.update(status='RUNNING',cut_authorized=False)
    (a.out/'status.json').write_text(json.dumps(manifest,indent=2))
    try:
        versions={t:subprocess.check_output([t,'--version'],text=True,stderr=subprocess.STDOUT).splitlines()[0]
                  for t in ('samtools','minimap2')}
        (a.out/'versions.json').write_text(json.dumps(versions,indent=2))
        if a.lane=='graph':
            result=graphs(a)
        else:
            job,folder,final,fa,regions=setup(a)
            manifest['assessment_sha256']=job['assessment_sha256']
            result=anchors(a,fa,regions) if a.lane=='anchors' else hifi(a,job,fa,regions) if a.lane=='hifi' else hic(a,job,folder,final,fa,regions)
        manifest.update(status='SUCCESS',result=result)
    except Exception as exc:
        manifest.update(status='FAILED',error=str(exc));raise
    finally:
        (a.out/'status.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':main()
