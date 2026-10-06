#!/usr/bin/env python3
"""Specific CTlk correction hypotheses; isolated copies, no production changes."""
import argparse
from collections import Counter, defaultdict
from bisect import bisect_left, bisect_right
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil

from run_ctlk_adjudication import CORES, paf, map_fasta, classify_anchor, peers
from trace_chimera_origins import Fasta
from chimera_origin import table, sha
from assess_junction_batch import write_rows
from junction_focus import intersection_length

TASKS = ['H01', 'J01', 'J02', 'J05', 'J06', 'J07']


def split_record(seq, start, end):
    if not 0 < start < end < len(seq) or set(seq[start:end].upper()) != {'N'}:
        raise ValueError('Requested split is not an internal all-N gap')
    # Preserve the entire gap as terminal Ns on the left record.
    left, right = seq[:end], seq[end:]
    if left+right != seq:
        raise ValueError('Sequence preservation failed')
    return left, right


def query_target_blocks(hit):
    """PAF cg coordinates including reverse query orientation and indels."""
    cg = next((t[5:] for t in hit['tags'].split(';') if t.startswith('cg:Z:')), None)
    if cg is None:
        raise ValueError('PAF missing CIGAR')
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cg)
    if ''.join(n+op for n, op in ops) != cg:
        raise ValueError('Invalid PAF CIGAR')
    t, q = hit['ts'], hit['qs'] if hit['strand'] == '+' else hit['qe']
    blocks = []
    for n, op in ops:
        n = int(n)
        if op in 'M=X':
            q2 = q+n if hit['strand'] == '+' else q-n
            blocks.append((t, t+n, min(q,q2), max(q,q2)))
        if op in 'MDN=X': t += n
        if op in 'MI=X': q += n if hit['strand'] == '+' else -n
    if t != hit['te'] or q != (hit['qe'] if hit['strand'] == '+' else hit['qs']):
        raise ValueError('PAF CIGAR and endpoints disagree')
    return blocks


def project_query(hit, position):
    for ts,te,qs,qe in query_target_blocks(hit):
        if qs <= position <= qe:
            return ts+position-qs if hit['strand']=='+' else te-(position-qs)
    return None


def source(a, hap, lane):
    folder=a.previous/f'{lane}_hap{hap}'
    status=json.loads((folder/'status.json').read_text())
    if status['status']!='SUCCESS':raise ValueError(f'Incomplete retained lane: {folder}')
    return folder,status


def assessment_fasta(a,hap):
    folder,status=source(a,hap,'anchors')
    origin=json.loads((folder/'origin_provenance.json').read_text())
    final=next(s for s in origin['stages'] if s['stage']=='assessment')
    if sha(Path(final['fasta']))!=status['assessment_sha256']:
        raise ValueError('Assessment FASTA checksum changed')
    return Fasta(Path(final['fasta']),a.out/'assessment.fa','samtools'),status


def h01(a):
    fa,status=assessment_fasta(a,1)
    start,end=63051225,63051325
    original=fa.fetch('scaffold_1',0,fa.lengths['scaffold_1'])
    left,right=split_record(original,start,end)
    digest=lambda s:hashlib.sha256(s.encode()).hexdigest()
    with (a.out/'unjoined.fa').open('w') as dest:
        for chrom,length in fa.lengths.items():
            if chrom=='scaffold_1':
                dest.write(f'>scaffold_1__left\n{left}\n>scaffold_1__right\n{right}\n')
            else:
                dest.write(f'>{chrom}\n{fa.fetch(chrom,0,length)}\n')
    import subprocess
    subprocess.run(['samtools','faidx',str(a.out/'unjoined.fa')],check=True)
    lengths={f[0]:int(f[1]) for f in (line.split('\t') for line in (a.out/'unjoined.fa.fai').read_text().splitlines())}
    expected={c:n for c,n in fa.lengths.items() if c!='scaffold_1'}
    expected.update(scaffold_1__left=len(left),scaffold_1__right=len(right))
    if lengths!=expected or sum(lengths.values())!=sum(fa.lengths.values()):
        raise ValueError('Written split assembly failed record/base accounting')
    write_rows(a.out/'coordinate_lift.tsv',[
        dict(old='scaffold_1',old_start=0,old_end=end,new='scaffold_1__left',new_start=0),
        dict(old='scaffold_1',old_start=end,old_end=len(original),new='scaffold_1__right',new_start=0)])
    folder,_=source(a,1,'anchors')
    # Relabel exact anchor coordinates; no peer remapping is required for identical sequence.
    annotations=[]
    with (folder/'anchor_correspondence.tsv').open() as src:
        for r in csv.DictReader(src,delimiter='\t'):
            if r['scaffold']!='scaffold_1':continue
            lo,hi=int(r['start']),int(r['end'])
            state='straddles_cut' if lo<end<hi else 'left' if hi<=end else 'right'
            annotations.append(dict(r,new_side=state,new_start=lo if state=='left' else lo-end if state=='right' else '.'))
    write_rows(a.out/'relabeled_anchors.tsv',annotations)
    # Reclassify existing endpoint bins. Boundary-straddling bins remain ambiguous.
    hicfolder,_=source(a,1,'hic'); counts=Counter()
    with gzip.open(hicfolder/'contact_bins.tsv.gz','rt') as src:
        for r in csv.DictReader(src,delimiter='\t'):
            if r['scaffold1']!='scaffold_1' and r['scaffold2']!='scaffold_1':continue
            sides=[]
            for suffix in ('1','2'):
                chrom=r['scaffold'+suffix];lo=int(r['bin'+suffix])*int(r['bin_bp']);hi=min(lo+int(r['bin_bp']),fa.lengths[chrom])
                sides.append(chrom if chrom!='scaffold_1' else 'left' if hi<=end else 'right' if lo>=end else 'boundary_bin')
            counts[r['library'],*sorted(sides)]+=int(r['count'])
    write_rows(a.out/'relabeled_contacts.tsv',[dict(library=l,side1=x,side2=y,count=n) for (l,x,y),n in sorted(counts.items())])
    for name in ('contact_scales.tsv','distance_profiles.tsv'):
        shutil.copy2(hicfolder/name,a.out/name)
    return dict(hypothesis='UNJOIN_UNSUPPORTED',assessment_sha256=status['assessment_sha256'],
        original_length=len(original),left_length=len(left),right_length=len(right),
        original_sequence_sha256=digest(original),reconstructed_sequence_sha256=digest(left+right),
        original_total_bases=sum(fa.lengths.values()),new_total_bases=sum(lengths.values()),
        original_records=len(fa.lengths),new_records=len(lengths),
        limitations='Isolated gap unjoin only; 100 Ns preserved at left terminus. Contact evidence unchanged and only relabeled; coarse boundary bins excluded from confident side assignment. No re-scaffolding or biological fusion verdict. Right piece still contains chr7/12 and requires separate adjudication.')


def graph_boundaries(a,hap,j,raw_start,raw_end):
    folder,_=source(a,hap,'graph'); endpoints=[]; coverage=[]
    inventory=table(folder/'graph_inventory.tsv')
    for treatment in ('both','hifi_only'):
        subset=[r for r in inventory if r['treatment']==treatment]
        # Runner numbers all native graphs per treatment in lexicographic order.
        for index,r in enumerate(subset):
            if not r['graph'].endswith('.r_utg.gfa') or r['status']!='mapped':continue
            dest=folder/f'{treatment}_{index}'
            for hit in table(dest/'seed_placements.tsv'):
                if hit['query']!=j:continue
                if int(hit['mapq'])<20 or int(hit['matches'])/int(hit['block'])<.999:continue
                offset=raw_start-100000
                for side,key in (('start','qs'),('end','qe')):
                    pos=offset+int(hit[key])
                    if raw_start<=pos<=raw_end:
                        endpoints.append(dict(treatment=treatment,segment=hit['target'],boundary_kind=side,
                                              raw_position=pos,query_start=hit['qs'],query_end=hit['qe']))
                typed={k:int(hit[k]) if k in ('qs','qe','ts','te') else hit[k] for k in hit}
                lo=project_query(typed,100000);hi=project_query(typed,int(hit['query_length'])-100000)
                if lo is None or hi is None:continue
                lo,hi=sorted((lo,hi));records=[]
                with (dest/'local.gfa').open() as src:
                    for line in src:
                        f=line.rstrip().split('\t')
                        if f[0]=='A' and f[1]==hit['target']:
                            records.append((int(f[2]),int(f[2])+int(f[6])-int(f[5]),f[4]))
                for pos in range(lo,hi,1000):
                    names={name for x,y,name in records if x<=pos<y}
                    coverage.append(dict(treatment=treatment,segment=hit['target'],segment_position=pos,
                                         graph_A_molecules_covering=len(names),interpretation='native offsets, not uniquely remapped read support'))
    write_rows(a.out/'raw_graph_endpoints.tsv',endpoints)
    write_rows(a.out/'native_graph_coverage.tsv',coverage)
    return endpoints


def evaluate_reads(a,hap,boundaries):
    folder,_=source(a,hap,'hifi')
    refs={r['name']:r for r in table(folder/'competitive_reference.tsv')}
    candidate=[]; controls=[]
    for r in table(folder/'regions.tsv'):
        if r['role']=='nearby_control':
            controls.append(dict(id=r['id'],scaffold=r['scaffold'],start=int(r['start']),end=int(r['end']),role='comparison_control'))
    tests=boundaries+controls
    indexed=defaultdict(list)
    for test in tests:indexed[test['scaffold']].append(test)
    indexes={}
    for chrom,rs in indexed.items():
        rs.sort(key=lambda r:r['start'])
        indexes[chrom]=(rs,[r['start'] for r in rs],max(r['end']-r['start'] for r in rs))
    # Preserve target/query coordinates for all competitor alignments of selected molecules.
    selected=set(); evidence=[]; placement_file=folder/'read_placements.paf'
    if not placement_file.is_file():raise ValueError('Retained competitive PAF missing')
    for hit in paf(placement_file):
        ref=refs[hit['target']]
        if f'Sde-CTlk_104_hap{hap}_round2_scaffolds.fa'!=ref['assembly']:continue
        blocks=query_target_blocks(hit)
        intervals=[(int(ref['start'])+x,int(ref['start'])+y) for x,y,_,_ in blocks]
        if not intervals or ref['scaffold'] not in indexes:continue
        rs,starts,width=indexes[ref['scaffold']]
        first=min(x for x,y in intervals);last=max(y for x,y in intervals)
        for test in rs[bisect_left(starts,first-2000-width):bisect_right(starts,last+2000)]:
            lo,hi=test['start'],test['end']
            left=intersection_length(intervals,[(lo-2000,lo)])
            right=intersection_length(intervals,[(hi,hi+2000)])
            if not left and not right:continue
            selected.add(hit['query'])
            evidence.append(dict(id=test['id'],role=test['role'],read=hit['query'],mapq=hit['mapq'],
                left_aligned_bp=left,right_aligned_bp=right,spans_1kb_each=left>=1000 and right>=1000,
                query_start=hit['qs'],query_end=hit['qe'],query_length=hit['query_length'],
                target_start=int(ref['start'])+hit['ts'],target_end=int(ref['start'])+hit['te'],
                matches=hit['matches'],block=hit['block'],tags=hit['tags']))
    write_rows(a.out/'narrow_read_evidence.tsv',evidence)
    with (a.out/'selected_molecule_competitors.tsv').open('w') as dest:
        fields=['read','assembly','scaffold','start','end','query_start','query_end','mapq','matches','block','tags']
        writer=csv.DictWriter(dest,fieldnames=fields,delimiter='\t');writer.writeheader()
        for hit in paf(placement_file):
            if hit['query'] not in selected:continue
            ref=refs[hit['target']]
            writer.writerow(dict(read=hit['query'],assembly=ref['assembly'],scaffold=ref['scaffold'],
                start=int(ref['start'])+hit['ts'],end=int(ref['start'])+hit['te'],query_start=hit['qs'],query_end=hit['qe'],
                mapq=hit['mapq'],matches=hit['matches'],block=hit['block'],tags=hit['tags']))
    summary=[]
    for test in tests:
        rs=[r for r in evidence if r['id']==test['id']]
        summary.append(dict(**test,spanning_molecules_any_mapq=len({r['read'] for r in rs if r['spans_1kb_each']}),
            spanning_molecules_mapq20=len({r['read'] for r in rs if r['spans_1kb_each'] and r['mapq']>=20}),
            placement_records=len(rs),interpretation='emitted competitive placements, ascertainment/clipping limits retained'))
    write_rows(a.out/'narrow_read_summary.tsv',summary)
    return summary


def boundary(a):
    hap,(_,chrom,lo,hi)=next((hap,core) for hap,cores in CORES.items() for core in cores if core[0]==a.task)
    fa,_=assessment_fasta(a,hap)
    raw=Fasta(a.assessment/'assembly/contig/hifiasm'/f'Sde-CTlk_104.hap{hap}.p_ctg.fasta',a.out/'raw.fa','samtools')
    offset=max(0,lo-250000);stop=min(raw.lengths[chrom],hi+250000)
    query=a.out/'raw_window.fa';query.write_text(f'>{a.task}\n{raw.fetch(chrom,offset,stop)}\n')
    hits=map_fasta(a,fa.path,query,a.out/'raw_to_assessment')
    exact=[h for h in hits if h['qs']==0 and h['qe']==stop-offset and h['matches']==h['block'] and h['mapq']>=20]
    if len(exact)!=1:raise ValueError('Raw window lacks one exact full assessment placement; no coordinate projection allowed')
    hit=exact[0]
    project=lambda rawpos:project_query(hit,rawpos-offset)
    endpoints=graph_boundaries(a,hap,a.task,lo,hi)
    proposals=[]
    if a.task=='J02':
        for rawpos in (4085098,4091627):
            pos=project(rawpos)
            if pos is None:raise ValueError('J02 seam projection failed')
            proposals.append(dict(id=f'J02_raw_{rawpos}',scaffold=hit['target'],start=pos,end=pos,role='graph_boundary_test'))
        projected=sorted((project(4085098),project(4091627)))
        proposals.append(dict(id='J02_seam_interval',scaffold=hit['target'],start=projected[0],end=projected[1],role='graph_seam_test'))
    for e in endpoints:
        pos=project(e['raw_position'])
        if pos is not None and all(r['start']!=pos for r in proposals):
            proposals.append(dict(id=f"{a.task}_raw_{e['raw_position']}",scaffold=hit['target'],start=pos,end=pos,role='graph_boundary_test'))
    # Dense smaller anchors supplement (not replace) the existing broad correspondences.
    center_lo,center_hi=sorted((project(lo),project(hi)))
    tilelo=max(0,center_lo-750000);tilehi=min(fa.lengths[hit['target']],center_hi+750000)
    meta=[]
    with (a.out/'tiles.fa').open('w') as dest:
        for size in (2000,5000,10000):
            for start in range(tilelo,tilehi-size+1,10000):
                seq=fa.fetch(hit['target'],start,start+size);q=f't{len(meta):05d}'
                meta.append(dict(query=q,scaffold=hit['target'],start=start,end=start+size,length=size,
                                 ambiguous_bases=sum(b not in 'ACGT' for b in seq.upper())))
                dest.write(f'>{q}\n{seq}\n')
    write_rows(a.out/'tile_coordinates.tsv',meta)
    correspondences=[]
    for i,peer in enumerate(peers(a)):
        print(f'{a.task}: smaller tiles against {peer.name}',flush=True)
        placements=map_fasta(a,peer,a.out/'tiles.fa',a.out/f'peer_{i}')
        grouped=defaultdict(list)
        for p in placements:grouped[p['query']].append(p)
        nm=a.assessment/'assembly/harmonization'/f"{peer.name.split('_round2')[0]}.harmonized_name_map.tsv"
        labels=defaultdict(set)
        for r in table(nm):labels[r['old_name']].add(r['new_name'])
        shutil.copy2(nm,a.out/f'peer_{i}.name_map.tsv')
        for m in meta:
            state,best=classify_anchor(grouped[m['query']],m['length'])
            if m['ambiguous_bases']:state='ambiguous_sequence'
            correspondences.append(dict(**m,peer=peer.name,status=state,target=best['target'] if best else '.',
                labels=';'.join(sorted(labels[best['target']])) if best else '.',
                target_start=best['ts'] if best else '.',target_end=best['te'] if best else '.',
                strand=best['strand'] if best else '.',mapq=best['mapq'] if best else '.'))
    write_rows(a.out/'smaller_anchor_correspondence.tsv',correspondences)
    # Always include the adjacent literal H01 gap in J01's comparison.
    if a.task=='J01':
        proposals.append(dict(id='H01_gap',scaffold='scaffold_1',start=63051225,end=63051325,role='adjacent_gap_test'))
    write_rows(a.out/'physical_boundary_tests.tsv',proposals)
    # Positions are assay probes, not cutting instructions or midpoint proposals.
    probes=[dict(id=f'{a.task}_probe_{pos}',scaffold=hit['target'],start=pos,end=pos,role='local_coverage_probe')
            for pos in range(center_lo,center_hi+1,1000)]
    write_rows(a.out/'local_coverage_probes.tsv',probes)
    evaluate_reads(a,hap,proposals+probes)
    graphfolder,_=source(a,hap,'graph')
    for name in ('graph_inventory.tsv',):shutil.copy2(graphfolder/name,a.out/name)
    return dict(hap=hap,raw_contig=chrom,raw_core=[lo,hi],assessment_core=[center_lo,center_hi],
        assessment_scaffold=hit['target'],raw_to_assessment_strand=hit['strand'],physical_boundary_tests=len(proposals),local_coverage_probes=len(probes),
        decision='REVIEW_REQUIRED',limitations='No cuts exported. Smaller anchors still require competing-placement checks and multiple independent tiles/individuals. Existing molecule placements are reused, not new full-read retrieval. Graph segment endpoints are hypotheses, not proven breakpoints. If no physical endpoints occur in the core, localization remains a correspondence question.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--task',choices=TASKS,required=True)
    p.add_argument('--previous',type=Path,required=True)
    p.add_argument('--assessment',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--threads',type=int,default=16)
    a=p.parse_args()
    if a.threads<=0:p.error('Positive threads required')
    a.out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(Path(__file__),a.out/'runner.py')
    status=dict(task=a.task,status='RUNNING',production_modified=False,cut_authorized=False)
    try:
        status['result']=h01(a) if a.task=='H01' else boundary(a)
        status['status']='SUCCESS'
    except Exception as exc:
        status.update(status='FAILED',error=str(exc));raise
    finally:
        (a.out/'status.json').write_text(json.dumps(status,indent=2)+'\n')


if __name__=='__main__':main()
