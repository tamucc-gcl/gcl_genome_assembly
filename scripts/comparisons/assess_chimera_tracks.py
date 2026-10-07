"""Offline assessment of chromosome tracks from retained PAFs; never authorizes cuts."""
import argparse,csv,html,io,json,re,tarfile,hashlib
from pathlib import Path
from collections import defaultdict,Counter


import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'py_scripts'))
from chimera_tracks import blocks,tracks,bins,side_summary


def assess(packet,out,mode="locus_unique"):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);summary=[];tile_rows=[];pages=[];all_data=[]
    colors=['#4477aa','#ee6677','#228833','#ccbb44','#66ccee','#aa3377','#bbbbbb','#332288','#88ccee','#44aa99','#117733','#999933','#ddcc77','#cc6677','#882255']
    palette={'chr'+str(n+1):c for n,c in enumerate(colors)}
    with tarfile.open(packet) as archive:
        members={m.name:m for m in archive.getmembers() if m.isfile()}
        for name in members:
            if not name.endswith('control_registry.json'):continue
            root=name.rsplit('/',1)[0]+'/';assembly=root.rstrip('/').split('/')[-1].replace('.sequence_context','')
            registry=json.load(archive.extractfile(members[name]));provenance=json.load(archive.extractfile(members[root+'provenance.json']))
            measurements=json.load(archive.extractfile(members[root+'decision_measurements.json']))
            gaps=list(csv.DictReader(io.StringIO(archive.extractfile(members[root+'literal_gap_catalog.tsv']).read().decode()),delimiter='\t'))
            candidates={k:v for k,v in registry['intervals'].items() if 'decision_id' in v};by_candidate=defaultdict(list)
            for index,peer in enumerate(provenance['peers']):
                grouped=defaultdict(list)
                for line in archive.extractfile(members[root+'peer_%d.paf'%index]).read().decode().splitlines():
                    f=line.split('\t')
                    if f[0] in candidates:grouped[f[0]].append(f)
                for key,interval in candidates.items():
                    track=tracks(grouped[key],peer.get('chromosome_labels',{}),mode);length=interval['end']-interval['start'];tiles=bins(track,length)
                    lo=interval['lo']-interval['start'];hi=interval['hi']-interval['start']
                    left=side_summary(track,tiles,0,lo);right=side_summary(track,tiles,hi,length)
                    qualified=left['qualified'] and right['qualified']
                    relation=('different_chromosomes' if left['chrom']!=right['chrom'] else 'same_chromosome') if qualified else 'uninformative'
                    band=None;nearby=[]
                    if relation=='different_chromosomes':
                        l=[b for b in tiles if b['hi']<=lo and b['chrom']==left['chrom']];r=[b for b in tiles if b['lo']>=hi and b['chrom']==right['chrom']]
                        band=[interval['start']+max(b['hi'] for b in l),interval['start']+min(b['lo'] for b in r)]
                        nearby=[g for g in gaps if g['scaffold']==interval['scaffold'] and band[0]<=int(g['lo'])<int(g['hi'])<=band[1]]
                    row=dict(assembly=assembly,candidate=key,peer=peer['id'],sample=peer['sample'],auto_evidence=peer['auto_evidence'],
                        left_chrom=left['chrom'],right_chrom=right['chrom'],left_aligned_bp=left['aligned_bp'],right_aligned_bp=right['aligned_bp'],
                        left_coverage=left['coverage'],right_coverage=right['coverage'],left_dominance=left['dominance'],right_dominance=right['dominance'],
                        left_bins=left['informative_bins'],right_bins=right['informative_bins'],relationship=relation,
                        band_start=band[0] if band else '',band_end=band[1] if band else '',gaps_in_band=len(nearby))
                    summary.append(row);by_candidate[key].append((row,tiles));all_data.append(dict(summary=row,tiles=tiles))
                    for b in tiles:tile_rows.append(dict(assembly=assembly,candidate=key,peer=peer['id'],scaffold=interval['scaffold'],start=interval['start']+b['lo'],end=interval['start']+b['hi'],chrom=b['chrom'] or '.',aligned_bp=b['aligned_bp'],ambiguous_bp=b['ambiguous_bp'],dominance=b['dominance']))
            for key,interval in candidates.items():
                rows=by_candidate[key];m=measurements[interval['decision_id']];length=interval['end']-interval['start'];x=lambda q:210+q/length*960
                title=assembly+' / '+('H01 location (discovered gap)' if m.get('gap_start')==63051225 else key)
                svg=['<svg viewBox="0 0 1200 330" role="img" aria-label="Chromosome evidence tracks">']
                for index,(row,tiles) in enumerate(rows):
                    y=35+index*27;svg.append('<text x="5" y="%d">%s</text>'%(y+12,html.escape(row['peer'])))
                    svg.append('<rect x="210" y="%d" width="960" height="18" fill="#f0f0f0"/>'%y)
                    for b in tiles:
                        if b['aligned_bp'] or b['ambiguous_bp']:
                            color=palette.get(b['chrom'],'#bcbcbc');opacity=max(.2,b['aligned_bp']/max(1,b['hi']-b['lo'])) if b['chrom'] else .6
                            svg.append('<rect x="%.2f" y="%d" width="%.2f" height="18" fill="%s" opacity="%.2f"><title>%s</title></rect>'%(x(b['lo']),y,(b['hi']-b['lo'])/length*960,color,opacity,html.escape(json.dumps(b))))
                for q in (interval['lo']-interval['start'],interval['hi']-interval['start']):svg.append('<line x1="%.2f" x2="%.2f" y1="20" y2="285" stroke="black" stroke-dasharray="4 3"/>'%(x(q),x(q)))
                svg.append('<text x="210" y="310">%.3f Mb</text><text x="1080" y="310">%.3f Mb</text></svg>'%(interval['start']/1e6,interval['end']/1e6))
                table='<table><tr><th>Peer</th><th>Left chromosome</th><th>Right chromosome</th><th>Left / right aligned kb</th><th>Dominance left / right</th><th>Assessment</th><th>Gaps in band</th></tr>'
                for r,b in rows:table+='<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [r['peer'],r['left_chrom'],r['right_chrom'],'%.1f / %.1f'%(r['left_aligned_bp']/1000,r['right_aligned_bp']/1000),'%.0f%% / %.0f%%'%(100*r['left_dominance'],100*r['right_dominance']),r['relationship'],r['gaps_in_band']])+'</tr>'
                table+='</table>'
                far='<table><tr><th>Offset kb</th><th>Library</th><th>Cross pairs</th><th>Within left/right</th><th>Matched sequence/gap controls</th><th>Candidate upper / control minimum</th><th>Informative?</th></tr>'
                for t in m.get('farther_contact_evidence',{}).get('trials',[]):
                    c=t['raw_counts'];p=t['control_populations'];ratio=t['upper_count_allowance_ratio']/t['minimum_control_ratio'] if t['minimum_control_ratio'] else None
                    far+='<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [t['offset_bp']//1000,t['library'],c.get('cross',0),str(c.get('left_within',0))+'/'+str(c.get('right_within',0)),str(p['continuous_control'])+'/'+str(p['gap_control']),round(ratio,3) if ratio is not None else 'unavailable',t['informative']])+'</tr>'
                far+='</table>'
                pages.append('<section><h2>'+html.escape(title)+'</h2><p>'+html.escape(interval['scaffold']+': '+str(interval['lo'])+'–'+str(interval['hi'])+' (0-based).')+'</p>'+''.join(svg)+table+'<h3>Existing farther-contact evidence (unchanged)</h3>'+far+'<p>Qualified immediate HiFi bridges: '+str(m.get('hifi_spanning_molecules'))+'. Absence over a broad interval is not seam evidence. Native graph contradiction: '+html.escape(str(m.get('graph_contradiction')))+'.</p></section>')
    for name,rows in [('peer_track_summary.tsv',summary),('chromosome_track_bins.tsv',tile_rows)]:
        with (out/name).open('w',newline='',encoding='utf-8') as handle:
            w=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (out/'track_assessment.json').write_text(json.dumps(dict(packet_sha256=hashlib.file_digest(open(packet,'rb'),'sha256').hexdigest(),diagnostic_only=True,assignment_mode=mode,records=all_data),indent=2),encoding='utf-8')
    legend=' '.join('<span style="color:%s">■ %s</span>'%(c,k) for k,c in palette.items())
    document='<!doctype html><meta charset="utf-8"><title>CTlk chromosome evidence assessment</title><style>body{font:15px system-ui;margin:30px;max-width:1400px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:6px;border:1px solid #ddd}svg{width:100%;font:12px system-ui}section{margin:40px 0;padding-top:15px;border-top:2px solid #aaa}</style><h1>Chromosome tracks: retrospective evidence assessment</h1><p>Diagnostic only. No automatic policy or assembly changes. Aggregates uniquely placed aligned bases across retained records; does not loosen MAPQ or identity filters. Overlaps count once. Competing near-equal placements are masked, including alternate loci on the same chromosome. Identity is record-wide because these archived CIGARs do not encode per-base mismatches.</p><p>50 kb display bins; color needs ≥10 kb assigned bases and ≥90% chromosome dominance. Opacity shows assigned coverage. Gray means ambiguous/mixed/insufficient support; pale background includes absent or filtered alignment. Dashed lines mark the interval edges, not validated cuts. Side qualification for this assessment needs ≥100 kb assigned bases, ≥90% dominance, and ≥3 bins. These are exploratory display criteria, not calibrated cut thresholds.</p><p>'+legend+'</p><p>Peer haplotypes are shown separately for inspection but are not independent individuals. Same-individual and poor-quality peers are contextual only. No alignment is not evidence of continuity or misassembly. Charts cover the retained local windows, not full chromosome lengths.</p>'+''.join(pages)
    document=document.replace('<h1>Chromosome tracks: retrospective evidence assessment</h1>','<h1>Chromosome tracks: '+html.escape(mode)+' assessment</h1>')
    if mode=='chromosome_unique':document=document.replace('Competing near-equal placements are masked, including alternate loci on the same chromosome.','Competing near-equal placements on different or unknown chromosomes are masked. Alternate loci confined to the same chromosome remain chromosome-assigned; their exact locus is unresolved. These assignments cannot demonstrate a molecular bridge across a join.')
    (out/'report.html').write_text(document,encoding='utf-8')
    print(json.dumps(dict(comparisons=len(summary),relationships=dict(Counter(r['relationship'] for r in summary)))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--packet',required=True);p.add_argument('--out',required=True);p.add_argument('--mode',choices=['locus_unique','chromosome_unique'],default='locus_unique');a=p.parse_args();assess(a.packet,a.out,a.mode)
