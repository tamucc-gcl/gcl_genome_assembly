#!/usr/bin/env python3
"""Focused read-only follow-up. Reuses HiFi placements; never edits assemblies."""
import argparse
from collections import Counter, defaultdict
import gzip
import itertools
import json
import re
from pathlib import Path
import shutil

from trace_chimera_origins import Fasta, map_queries, run
from chimera_origin import table, sha, agp_rows
from assess_junction_batch import write_rows, pairs_dictionary, verify_agp
from junction_assessment import Contacts, Lift, library_for
from junction_focus import read_blocks, focal_query, intersection_length, anchor_intervals, anchor_pair


def alternatives(packet, dest, folder, region, names):
    """Two passes over retained compact SAM placements; no new read mapping."""
    mapped = names[region['assembly'],region['scaffold']]
    path = packet/folder/'competitive_alignments.tsv.gz'
    if not path.is_file():
        raise ValueError('Retained competitive placements missing')
    chosen = {}
    def records():
        with gzip.open(path,'rt') as src:
            next(src)
            for line in src:
                f=line.rstrip('\n').split('\t')
                yield dict(qname=f[0],flag=int(f[1]),rname=f[2],pos=int(f[3])-1,
                           mapq=int(f[4]),cigar=f[5],tags=dict(t.split(':',2)[::2] for t in f[6].split(';') if t.count(':')>=2))
    for r in records():
        if r['flag'] & 0xF04 or r['rname'] != mapped:
            continue
        length, blocks=read_blocks(r['pos'],r['cigar'],bool(r['flag']&16))
        core=focal_query(blocks,int(region['start']),int(region['end']),bool(r['flag']&16))
        bases=sum(b-a for a,b in core)
        if bases >= 1000:
            if r['qname'] in chosen:
                raise ValueError('Multiple retained primary records for one molecule')
            chosen[r['qname']]=(r,length,core,bases)
    output=[]
    for r in records():
        if r['qname'] not in chosen or r['flag']&4 or r['cigar']=='*':
            continue
        primary,length,core,bases=chosen[r['qname']]
        altlength,blocks=read_blocks(r['pos'],r['cigar'],bool(r['flag']&16))
        if altlength != length:
            raise ValueError('Read length differs across retained placements')
        fraction=intersection_length(core,[(c,d) for _,_,c,d in blocks])/bases
        other=next((k for k,v in names.items() if v==r['rname']),('.',r['rname']))
        output.append(dict(region=region['id'],read=r['qname'],primary_mapq=primary['mapq'],
            primary_AS=primary['tags'].get('AS','.'),primary_focal_aligned_bp=bases,
            assembly=other[0],scaffold=other[1],start=r['pos'],end=max((b for _,b,_,_ in blocks),default=r['pos']),
            flag=r['flag'],mapq=r['mapq'],AS=r['tags'].get('AS','.'),NM=r['tags'].get('NM','.'),
            focal_read_fraction_aligned=fraction,cigar=r['cigar'],
            relationship='primary' if not r['flag']&0x900 else
                ('other_haplotype' if other[0]!=region['assembly'] else 'same_assembly_alternative')))
    write_rows(dest/f"{region['id']}.read_alternatives.tsv",output)
    per_read=[]
    grouped=defaultdict(list)
    for r in output:
        if r['relationship']!='primary' and r['focal_read_fraction_aligned']>=0.8:
            grouped[r['read']].append(r)
    for name,(primary,_,_,bases) in chosen.items():
        rs=grouped[name]
        per_read.append(dict(region=region['id'],read=name,primary_mapq=primary['mapq'],focal_aligned_bp=bases,
            emitted_alternatives_covering_80pct=len(rs),other_haplotype=sum(r['relationship']=='other_haplotype' for r in rs),
            same_assembly_alternatives=sum(r['relationship']=='same_assembly_alternative' for r in rs)))
    write_rows(dest/f"{region['id']}.read_summary.tsv",per_read)
    return dict(id=region['id'],source_placements_sha256=sha(path),primary_reads_selected=len(chosen),
                reads_with_emitted_alternative=sum(bool(grouped[n]) for n in chosen),
                limitation='Only originally selected reads and emitted placements; 80% overlap is descriptive, not a uniqueness or phase test. AS is whole-alignment score.')


def contacts(a, job, gaps, out):
    meta=job['hic']
    if meta.get('status')!='measured':
        raise ValueError('Source packet lacks measured Hi-C provenance')
    folder=a.packet/job['folder']
    origin=json.loads((folder/'origin_provenance.json').read_text())
    final=next(s for s in origin['stages'] if s['stage']=='assessment')
    if sha(Path(final['agp']))!=meta['agp_sha256'] or sha(Path(meta['source_fasta']))!=meta['source_fasta_sha256']:
        raise ValueError('Retained Hi-C source/AGP changed')
    agp=agp_rows(final['agp'])
    verify_agp(agp,job['fasta'])
    source_lengths={f[0]:int(f[1]) for f in (l.split('\t') for l in (folder/'hic_input.fa.fai').read_text().splitlines())}
    if pairs_dictionary(Path(meta['pairs']))!=source_lengths:
        raise ValueError('Pairs dictionary changed')
    manifest=folder/'hic_readsets.tsv'
    if sha(manifest)!=meta['readsets_sha256']:
        raise ValueError('Read-set manifest changed')
    libraries=table(manifest)
    lift=Lift(agp)
    counters={size:Contacts(gaps,job['fasta'].lengths,size) for size in a.contact_sizes}
    input_count=0
    with gzip.open(meta['pairs'],'rt') as src:
        for line in src:
            if line.startswith('#') or not line.strip():continue
            f=line.rstrip().split('\t')
            if len(f)<8:raise ValueError('Malformed pairs')
            input_count+=1
            if f[7]!='UU':continue
            x,y=lift.locate(f[1],int(f[2])-1),lift.locate(f[3],int(f[4])-1)
            if x is None or y is None:continue
            lib=library_for(f[0],libraries)
            for counter in counters.values():counter.add(lib,x,y)
    if input_count!=meta['audits']['input_pairs']:
        raise ValueError('Retained pairs count changed')
    rows,partners=[],[]
    for size,c in counters.items():
        sequence_windows={}
        for g in gaps:
            for side in ('left','right'):
                chrom,lo,hi=c.windows[g['id'],side]
                seq=job['fasta'].fetch(chrom,lo,hi) if lo<hi else ''
                sequence_windows[g['id'],side]=(sum(b not in 'ACGT' for b in seq),
                    sum(h['id']!=g['id'] and h['scaffold']==chrom and h['start']<hi and lo<h['end'] for h in job['all_gaps']))
        for lib in sorted({r['library_id'] for r in libraries}|set(c.total)):
            for g in gaps:
                rows.append(dict(id=g['id'],role=g['role'],flank_bp=size,library=lib,
                    crossing_pairs=c.cross[lib,g['id']],projected_UU_pairs=c.total[lib],
                    left_window_bp=c.windows[g['id'],'left'][2]-c.windows[g['id'],'left'][1],
                    right_window_bp=c.windows[g['id'],'right'][2]-c.windows[g['id'],'right'][1],
                    left_ambiguous_bases=sequence_windows[g['id'],'left'][0],right_ambiguous_bases=sequence_windows[g['id'],'right'][0],
                    left_other_gaps=sequence_windows[g['id'],'left'][1],right_other_gaps=sequence_windows[g['id'],'right'][1],
                    left_external_contacts=c.margin[lib,g['id'],'left'],right_external_contacts=c.margin[lib,g['id'],'right']))
                for side in ('left','right'):
                    hits=[(k,v) for k,v in c.partners.items() if k[:3]==(lib,g['id'],side)]
                    for k,n in sorted(hits,key=lambda kv:(-kv[1],kv[0]))[:20]:
                        partners.append(dict(id=g['id'],flank_bp=size,library=lib,side=side,
                            partner_scaffold=k[3],partner_start=k[4]*size,
                            partner_end=min((k[4]+1)*size,job['fasta'].lengths[k[3]]),count=n))
    write_rows(out/'contact_scales.tsv',rows)
    write_rows(out/'alternative_contact_bins.tsv',partners)
    # Reproduce the original 100 kb results exactly before interpreting new scales.
    original=table(folder/'hic_by_library.tsv')
    lookup={(r['id'],r['library']):r for r in original}
    for r in rows:
        if r['flank_bp']==100000:
            old=lookup[r['id'],r['library']]
            for key in ('crossing_pairs','projected_UU_pairs','left_external_contacts','right_external_contacts'):
                if int(old[key])!=r[key]:raise ValueError('100 kb replay differs from source packet')


def anchors(a,jobs):
    dest=a.out/'anchors';dest.mkdir()
    records=[]
    with (dest/'anchors.fa').open('w') as fasta:
        for job in jobs:
            for g in job['gaps']:
                for offset in a.offsets:
                    for side,iv in anchor_intervals(g['start'],g['end'],job['fasta'].lengths[g['scaffold']],a.anchor_size,offset).items():
                        ident=f"{g['id']}_{side}_{offset}"
                        r=dict(id=ident,gap=g['id'],assembly=job['assembly'],scaffold=g['scaffold'],side=side,offset=offset,status='outside_scaffold',start='.',end='.',ambiguous_bases='.',other_gaps_to_join='.')
                        if iv:
                            seq=job['fasta'].fetch(g['scaffold'],*iv)
                            lo,hi=(iv[0],g['start']) if side=='left' else (g['end'],iv[1])
                            n=sum(h['id']!=g['id'] and h['scaffold']==g['scaffold'] and h['start']<hi and lo<h['end'] for h in job['all_gaps'])
                            r.update(start=iv[0],end=iv[1],status='extracted',ambiguous_bases=sum(b not in 'ACGT' for b in seq),other_gaps_to_join=n)
                            fasta.write(f'>{ident}\n{seq}\n')
                        records.append(r)
    write_rows(dest/'windows.tsv',records)
    results=[];comparisons=[]
    for i,s in enumerate(table(a.packet/'sequence_comparisons/sources.tsv')):
        path=Path(s['path'])
        if sha(path)!=s['sha256']:raise ValueError('Peer FASTA changed')
        print(f"Mapping flank anchors to {path.name}",flush=True)
        target=Fasta(path,dest/f'peer_{i}.fa',a.samtools)
        hits=map_queries(a,target,dest/'anchors.fa',dest/f'peer_{i}_anchors')
        qualified=defaultdict(list)
        for h in hits:
            aligned=sum(int(n) for n,op in re.findall(r'(\d+)([=XMID])',h['cigar']) if op in '=XM')
            frac=aligned/h['query_length']
            results.append(dict(peer=i,peer_path=s['path'],**h,aligned_query_fraction=frac))
            if frac>=0.8:qualified[h['query']].append(h)
        for job in jobs:
            for g in job['gaps']:
                for offset in a.offsets:
                    left=qualified[f"{g['id']}_left_{offset}"]
                    right=qualified[f"{g['id']}_right_{offset}"]
                    for hs in (left,right):hs.sort(key=lambda h:(-h['mapq'],-h['matches'],h['target'],h['target_start']))
                    base=dict(gap=g['id'],peer=i,peer_path=s['path'],offset=offset,left_hits=len(left),right_hits=len(right),
                              combinations_truncated=len(left)>5 or len(right)>5)
                    if not left or not right:
                        comparisons.append(dict(base,status='missing_qualifying_flank'))
                    for l,r in itertools.product(left[:5],right[:5]):
                        status,distance=anchor_pair(l,r)
                        comparisons.append(dict(base,status=status,left_target=l['target'],right_target=r['target'],
                            left_mapq=l['mapq'],right_mapq=r['mapq'],left_strand=l['strand'],right_strand=r['strand'],
                            target_separation=distance,expected_assembly_separation=g['end']-g['start']+2*offset))
    write_rows(dest/'placements.tsv',results)
    fields=list(dict.fromkeys(k for row in comparisons for k in row))
    write_rows(dest/'flank_correspondence.tsv',comparisons,fields)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--packet',type=Path,required=True)
    p.add_argument('--focus',type=Path,required=True,help='TSV: id, assembly, scaffold, start, end, kind, assessment_sha256. kind is gap or read_region')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--threads',type=int,default=8)
    p.add_argument('--index-bases',type=int,default=8000000000)
    p.add_argument('--anchor-size',type=int,default=10000)
    p.add_argument('--offsets',type=int,nargs='+',default=[0,25000,100000,250000])
    p.add_argument('--contact-sizes',type=int,nargs='+',default=[25000,100000,250000,500000])
    p.add_argument('--samtools',default='samtools');p.add_argument('--minimap2',default='minimap2')
    a=p.parse_args()
    if min(a.threads,a.index_bases,a.anchor_size,*a.contact_sizes)<=0 or min(a.offsets)<0 or 100000 not in a.contact_sizes:
        p.error('Invalid settings; include 100000 in --contact-sizes for replay validation')
    focus=table(a.focus)
    from junction_assessment import validate_registry
    validate_registry([dict(r,sample='registry') for r in focus])
    if not focus or any(r['kind'] not in ('gap','read_region') for r in focus):p.error('Empty/invalid focus registry')
    status=json.loads((a.packet/'status.json').read_text())
    if status['status']!='SUCCESS':raise ValueError('Source packet did not complete')
    a.out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(a.focus,a.out/'focus.tsv');shutil.copy2(a.packet/'status.json',a.out/'source_status.json')
    shutil.copy2(a.packet/'sequence_comparisons/sources.tsv',a.out/'peer_sources.tsv')
    for src in (Path(__file__),Path(__file__).parents[2]/'py_scripts/junction_focus.py'):
        shutil.copy2(src,a.out/src.name)
    (a.out/'parameters.json').write_text(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},indent=2))
    for tool in (a.samtools,a.minimap2):run([tool,'--version'],stdout=a.out/(Path(tool).name+'.version.txt'))
    jobs=[];reads=[]
    names={(r['assembly'],r['scaffold']):r['mapped_name'] for r in table(a.packet/'competitive_reference_names.tsv')}
    for job in status['assemblies']:
        selected=[r for r in focus if r['assembly']==job['assembly']]
        if not selected:continue
        folder=a.packet/job['folder'];dest=a.out/job['folder'];dest.mkdir()
        source=json.loads((folder/'origin_provenance.json').read_text())
        final=next(s for s in source['stages'] if s['stage']=='assessment')
        digest=sha(Path(final['fasta']))
        if digest!=job['assessment_sha256'] or any(r['assessment_sha256']!=digest for r in selected):raise ValueError('Assessment checksum mismatch')
        job['fasta']=Fasta(Path(final['fasta']),dest/'assessment.fa',a.samtools)
        gaps=table(folder/'selected_gaps.tsv')
        for g in gaps:g.update(start=int(g['start']),end=int(g['end']))
        job['all_gaps']=table(folder/'all_verified_gaps.tsv')
        for g in job['all_gaps']:g.update(start=int(g['start']),end=int(g['end']))
        by_id={g['id']:g for g in gaps}
        for r in selected:
            if r['scaffold'] not in job['fasta'].lengths or int(r['end'])>job['fasta'].lengths[r['scaffold']]:raise ValueError('Focus outside assessment')
            if r['kind']=='gap':
                g=by_id.get(r['id'])
                if not g or any(str(g[k])!=r[k] for k in ('scaffold','start','end')):raise ValueError('Focus gap differs from verified packet')
            else:reads.append(alternatives(a.packet,dest,job['folder'],r,names))
        wanted={r['id'] for r in selected if r['kind']=='gap'}
        # Preserve all old comparison gaps on the same scaffold, not just the best-looking controls.
        scaffolds={by_id[x]['scaffold'] for x in wanted}
        job['gaps']=[g for g in gaps if g['id'] in wanted or (g['role']=='nearby_control' and g['scaffold'] in scaffolds)]
        if job['gaps']:
            print(f"Streaming retained Hi-C for {job['assembly']}",flush=True)
            write_rows(dest/'gaps.tsv',job['gaps']);contacts(a,job,job['gaps'],dest)
            jobs.append(job)
    if set(r['assembly'] for r in focus)-set(j['assembly'] for j in status['assemblies']):raise ValueError('Unknown assembly')
    anchors(a,jobs)
    (a.out/'status.json').write_text(json.dumps(dict(status='SUCCESS',cut_authorized=False,read_regions=reads,
        limitations='Descriptive comparisons; no calibrated misjoin likelihood. Anchors may cross other gaps: inspect ambiguous_bases. Self/historical alignments are not independent biological validation.'),indent=2))


if __name__=='__main__':main()
