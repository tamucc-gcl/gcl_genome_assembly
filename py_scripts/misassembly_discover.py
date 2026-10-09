#!/usr/bin/env python3
"""Construct cohort labels and discover chromosome-scale transitions in every assembly."""
import argparse
import json
from collections import defaultdict
from pathlib import Path
from misassembly_core import label_catalog,cohort_labels,independent_labels,read_paf,discover,sha,POLICY,coarse_tracks
from chimera_tracks import tracks


def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    c=sub.add_parser('catalog');c.add_argument('--manifest',required=True);c.add_argument('--alignments',required=True);c.add_argument('--out',required=True)
    d=sub.add_parser('discover')
    for name in ('catalog','assembly','fasta','alignments','out'):d.add_argument('--'+name,required=True)
    d.add_argument('--minimum-arm-bp',type=int,default=1000000)
    a=p.parse_args()
    if a.command=='catalog':
        records=json.loads(a.manifest)
        Path(a.out).write_text(json.dumps(cohort_labels(label_catalog(records),json.loads(a.alignments))),encoding='utf-8');return
    catalog=json.loads(Path(a.catalog).read_text());own=catalog[a.assembly];scope=set(own['scope'])
    peer_tracks=[];alignment_provenance=[];coverage=[]
    for alignment in json.loads(a.alignments):
        if a.assembly not in (alignment['a'],alignment['b']):raise ValueError('Unrelated pair alignment staged')
        reverse=a.assembly==alignment['b'];peer=catalog[alignment['a'] if reverse else alignment['b']]
        if peer['taxid']!=own['taxid']:raise ValueError('Cross-species chimera comparison')
        grouped=defaultdict(list)
        for f in read_paf(alignment['path'],reverse):
            if f[0] in scope:grouped[f[0]].append(f)
        labels=independent_labels(peer,own['sample'])
        mapped={sc:tracks(hits,labels,mode='chromosome',min_mapq=20,minimum_target=0) for sc,hits in grouped.items()}
        peer_tracks.append((peer,mapped))
        for scaffold in scope:
            assigned=sum(b['hi']-b['lo'] for b in mapped.get(scaffold,[]) if b.get('chrom'))
            length=own['scope_lengths'].get(scaffold) or (int(grouped[scaffold][0][1]) if grouped.get(scaffold) else None)
            coverage.append(dict(peer=peer['id'],sample=peer['sample'],eligible=peer['eligible'],scaffold=scaffold,assigned_bp=assigned,assigned_fraction=assigned/length if length else None))
        alignment_provenance.append(dict(a=alignment['a'],b=alignment['b'],sha256=sha(alignment['path'])))
    events,minor=discover(peer_tracks,own['sample'],a.minimum_arm_bp)
    for event in events:
        event['peer_tracks']=[dict(peer=peer['id'],sample=peer['sample'],eligible=peer['eligible'],
                blocks=[dict(b,lo=max(0,event['lo']-1000000,b['lo']),hi=min(event['hi']+1000000,b['hi'])) for b in mapped.get(event['scaffold'],[])
                        if b['lo']<event['hi']+1000000 and b['hi']>event['lo']-1000000]) for peer,mapped in peer_tracks]
    status='assessed' if scope and peer_tracks else 'no_chromosome_scope' if not scope else 'no_peer_alignments'
    informative=[peer['id'] for peer,mapped in peer_tracks if peer['eligible'] and peer['sample']!=own['sample'] and
                 any(sum(b['hi']-b['lo'] for b in track if b.get('chrom'))>=a.minimum_arm_bp for track in mapped.values())]
    if status=='assessed' and not informative:status='no_informative_peers'
    out=Path(a.out);out.mkdir()
    result=dict(policy=POLICY,assembly=a.assembly,sample=own['sample'],taxid=own['taxid'],assessment_sha256=sha(a.fasta),
                coordinate_stage='pre_finishing',status=status,informative_peers=informative,scope=sorted(scope),minimum_arm_bp=a.minimum_arm_bp,scope_unresolved=own['scope_unresolved'],events=events,minor_signals=minor,alignments=alignment_provenance,
                screening_coverage=coverage,labels={peer['id']:dict(labelled_scaffolds=sum(bool(v) for v in peer['chromosome_labels'].values()),composite_blocks=sum(len(v) for v in peer['chromosome_labels'].values() if isinstance(v,list))) for peer,_ in peer_tracks})
    result['coarse_context']=[dict(peer=peer['id'],sample=peer['sample'],eligible=peer['eligible'],scaffolds={sc:coarse_tracks(track) for sc,track in mapped.items()}) for peer,mapped in peer_tracks]
    (out/'discovery.json').write_text(json.dumps(result),encoding='utf-8')


if __name__=='__main__':main()
