"""Create human evidence packets and editable, unselected review rows. No auto decisions."""
import argparse,csv,hashlib,html,json
from pathlib import Path
from chimera_markdown import render, decision_context, assess_transitions
FIELDS=['selected','id','assembly','coordinate_stage','assessment_sha256','scaffold','action','cut_bp','gap_start','gap_end','review_start','review_end','review_range','localization_explanation',
        'decision_source','evidence_packet_id','reviewer','reason','localization_status','evidence_summary','report_path','source_candidate_id','chromosome_context','evidence_for_cut','evidence_against_cut','review_priority','evidence_limits','detected_transition','assessment_status','transition_id','preferred_candidate','bridge_status','related_candidate']

def write_table(path,rows):
    with path.open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=FIELDS,delimiter='\t');w.writeheader();w.writerows(rows)

def generate(assembly,calls,context,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);context=Path(context)
    provenance=json.loads((context/'provenance.json').read_text()) if (context/'provenance.json').exists() else {}
    measurements=json.loads((context/'decision_measurements.json').read_text()) if (context/'decision_measurements.json').exists() else {}
    with open(calls,encoding='utf-8') as h:records=list(csv.DictReader((l for l in h if not l.startswith('#')),delimiter='\t'))
    candidates={hashlib.sha256(json.dumps(r,sort_keys=True).encode()).hexdigest()[:20]:r for r in records}
    for key,m in measurements.items():
        if m.get('review_only'):candidates[key]=dict(scaffold=m['scaffold'],assembly_sha256=m['assessment_sha256'])
    rows=[];sections=[];flat=[]
    palette={('chr'+str(i+1)):color for i,color in enumerate(['#4477aa','#ee6677','#228833','#ccbb44','#66ccee','#aa3377','#bbbbbb','#332288','#88ccee','#44aa99','#117733','#999933','#ddcc77','#cc6677','#882255'])}
    def location(item):
        key,c=item;m=measurements.get(key,{})
        interval=provenance.get('intervals',{}).get(m.get('packet_interval_id'),{})
        scaffold=c['scaffold']
        digits=''.join(x for x in scaffold if x.isdigit())
        return (int(digits) if digits else 0,scaffold,int(m.get('cut_bp') or c.get('transition_lo') or interval.get('lo') or 0),key)
    for number,(key,c) in enumerate(sorted(candidates.items(),key=location),1):
        m=measurements.get(key,{});gap=m.get('verified_gap') is True
        coordinate=m.get('cut_bp','') if gap else ''
        summary='Literal gap proposed for human review' if gap else 'Transition unlocalized; exact internal cut requires explicit review'
        display_id='C%02d'%number
        row=dict.fromkeys(FIELDS,'');row.update(source_candidate_id=key,selected='NO',id=display_id,assembly=assembly,coordinate_stage='pre_finishing',
            assessment_sha256=c.get('assembly_sha256',m.get('assessment_sha256','')),scaffold=c['scaffold'],
            action='UNJOIN_UNSUPPORTED' if gap else 'UNRESOLVED',cut_bp=coordinate,gap_start=m.get('gap_start',''),gap_end=m.get('gap_end',''),
            decision_source='review',evidence_packet_id=assembly+'.sequence_context',localization_status='proposed_gap' if gap else 'unlocalized',
            evidence_summary=summary,report_path=assembly+'.review/report.md#candidate-'+display_id.lower())
        row.update(decision_context(m,provenance.get('sample',assembly.rsplit('_hap',1)[0])))
        interval=provenance.get('intervals',{}).get(m.get('packet_interval_id'),{})
        lo=m.get('gap_start') if gap else c.get('transition_lo',interval.get('lo'))
        hi=m.get('gap_end') if gap else c.get('transition_hi',interval.get('hi'))
        row.update(review_start=lo if lo is not None else '',review_end=hi if hi is not None else '',
                   review_range=('%s–%s'%(lo,hi)) if lo is not None and hi is not None else 'Interval unavailable')
        row['localization_explanation']=('Verified all-N gap; an exact gap-end cut is available' if gap else
            'Local HiFi continuity supports the sampled transition; chromosome assignments do not identify a failed seam' if m.get('local_path_support')=='supported_grid' else
            'Repeat or ambiguous alignment obscures the seam; the chromosome-transition interval remains available for review' if m.get('repeat_obscured_localization') else
            'Chromosome-transition interval is measured, but no unique failed seam or verified gap has been established inside it')
        rows.append(row);escaped=html.escape
        section='<section id="'+key+'"><h2>'+escaped(c['scaffold']+' / '+key)+'</h2><p>'+summary+'. Selected: NO. Proposed cut: '+str(coordinate or 'none')+'. All coordinates are 0-based on the original assessed FASTA.</p>'
        tracks=m.get('chromosome_tracks',[]);all_bins=[b for t in tracks for b in t['bins']];length=max((b['hi'] for b in all_bins),default=1)
        svg=['<svg viewBox="0 0 1200 %d" role="img" aria-label="Peer chromosome tracks">'%(max(80,40+len(tracks)*28))]
        for i,t in enumerate(tracks):
            y=20+i*28;role='same individual' if t['sample']==assembly.rsplit('_hap',1)[0] else 'independent eligible' if t['auto_evidence'] else 'context only'
            svg.append('<text x="0" y="%d">%s</text>'%(y+13,escaped(t['peer'])))
            svg.append('<rect x="210" y="%d" width="980" height="18" fill="#f3f3f3"/>'%y)
            for b in t['bins']:
                if b['aligned_bp'] or b['ambiguous_bp']:
                    color=palette.get(b['chrom'],'#aaa');opacity=max(.2,b['aligned_bp']/max(1,b['hi']-b['lo'])) if b['chrom'] else .6
                    svg.append('<rect x="%.2f" y="%d" width="%.2f" height="18" fill="%s" opacity="%.2f"><title>%s</title></rect>'%(210+b['lo']/length*980,y,(b['hi']-b['lo'])/length*980,color,opacity,escaped(json.dumps(b))))
            flat.append(dict(candidate=key,peer=t['peer'],sample=t['sample'],role=role,relationship=t['relationship'],
                left_chrom=t['left']['chrom'] or '.',right_chrom=t['right']['chrom'] or '.',left_aligned_bp=t['left']['aligned_bp'],right_aligned_bp=t['right']['aligned_bp'],
                left_coverage=t['left']['coverage'],right_coverage=t['right']['coverage'],left_dominance=t['left']['dominance'],right_dominance=t['right']['dominance']))
        section+=''.join(svg)+'</svg><p>Colors show chromosome assignment; opacity shows assigned coverage. Gray indicates ambiguity, mixed identity, or insufficient evidence. Pale background includes absent or filtered alignment. Local windows only; tracks do not prove molecular continuity. Different chromosomes means the assessed left and right sequences map to separate chromosomes in the peer, not that the peer has a join.</p>'
        interval=provenance.get('intervals',{}).get(m.get('packet_interval_id'),{})
        if interval:
            # Add original assessment positions; dashed edges describe hypotheses, not approval.
            lines=''.join('<line x1="%.2f" x2="%.2f" y1="5" y2="%d" stroke="black" stroke-dasharray="4 3"/>'%(210+(q-interval['start'])/length*980,210+(q-interval['start'])/length*980,20+len(tracks)*28) for q in [interval['lo'],interval['hi']])
            section=section.replace('</svg>',lines+'</svg><p>Displayed assessment window: '+str(interval['start'])+'–'+str(interval['end'])+' bp. Dashed lines mark hypothesis edges.</p>')
        by_individual={}
        for t in tracks:
            if t['auto_evidence'] and t['sample']!=assembly.rsplit('_hap',1)[0]:
                by_individual.setdefault(t['sample'],set())
                if t['relationship']!='uninformative':by_individual[t['sample']].add((t['relationship'],t['left']['chrom'],t['right']['chrom']))
        different=[sample for sample,pairs in by_individual.items() if len(pairs)==1 and next(iter(pairs))[0]=='different_chromosomes']
        same=[sample for sample,pairs in by_individual.items() if len(pairs)==1 and next(iter(pairs))[0]=='same_chromosome']
        missing=[sample for sample,pairs in by_individual.items() if not pairs]
        conflicting=[sample for sample,pairs in by_individual.items() if len(pairs)>1]
        row['evidence_summary']=row['review_priority']+'; '+row['chromosome_context']
        section+='<p><b>Independent chromosome context:</b> different chromosomes in '+escaped(', '.join(different) or 'none')+'; same chromosome in '+escaped(', '.join(same) or 'none')+'; missing in '+escaped(', '.join(missing) or 'none')+'; conflicting in '+escaped(', '.join(conflicting) or 'none')+'. Missing measurements neither support nor oppose a cut.</p>'
        table='<table><tr><th>Peer / role</th><th>Assessed left / right chromosome assignments</th><th>Aligned kb left/right</th><th>Chromosome context</th></tr>'
        for t,r in zip(tracks,[r for r in flat if r['candidate']==key]):
            table+='<tr>'+''.join('<td>'+escaped(str(v))+'</td>' for v in [t['peer']+' / '+r['role'],r['left_chrom']+' / '+r['right_chrom'],'%.1f / %.1f'%(r['left_aligned_bp']/1000,r['right_aligned_bp']/1000),t['relationship']])+'</tr>'
        section+=table+'</table>'
        link='../../sequence_context/'+assembly+'.sequence_context/'
        packet=m.get('packet_interval_id')
        if packet:
            for suffix,title in [('.controls.png','Immediate measurements and controls'),('.farther_contacts.png','Farther contact measurements')]:
                if (context/(packet+suffix)).exists():section+='<h3>'+title+'</h3><img src="'+escaped(link+packet+suffix)+'" alt="'+title+'">'
            section+='<p><a href="'+escaped(link+packet+'.igv.xml')+'">IGV: original assessment coordinates</a></p>'
        section+='<p>Qualified immediate HiFi spanning molecules: '+str(m.get('hifi_spanning_molecules','unavailable'))+'. Zero over a broad interval is not a failed seam. Graph context: '+escaped(str(m.get('graph_status','unavailable')))+'.</p>'
        hifi=m.get('hifi_raw',{})
        section+='<p>Immediate qualified HiFi flank molecules, left / right: '+str(hifi.get('left_molecules','unavailable'))+' / '+str(hifi.get('right_molecules','unavailable'))+'. Poor flank observability weakens absence-of-support evidence.</p>'
        contacts='<table><tr><th>Offset kb</th><th>Library</th><th>Cross pairs</th><th>Within left / right</th><th>Matched sequence / gap controls</th><th>Assay informative?</th></tr>'
        for trial in m.get('farther_contact_evidence',{}).get('trials',[]):
            raw=trial['raw_counts'];pop=trial['control_populations']
            contacts+='<tr>'+''.join('<td>'+escaped(str(v))+'</td>' for v in [trial['offset_bp']//1000,trial['library'],raw.get('cross',0),str(raw.get('left_within',0))+' / '+str(raw.get('right_within',0)),str(pop['continuous_control'])+' / '+str(pop['gap_control']),trial['informative']])+'</tr>'
        section+='<h3>Per-library contact measurements</h3>'+contacts+'</table><p>Control thresholds describe assay calibration, not permission to cut. Scaffolding Hi-C reads are corroboration, not an independent validation dataset.</p>'
        section+='<details><summary>Measured support, opposing and missing evidence</summary><pre>'+escaped(json.dumps(m,indent=2))+'</pre></details></section>'
        sections.append(section)
    assess_transitions(rows,measurements,provenance,candidates)
    write_table(out/'review.tsv',rows)
    if flat:
        with (out/'peer_track_summary.tsv').open('w',newline='',encoding='utf-8') as h:
            w=csv.DictWriter(h,fieldnames=list(flat[0]),delimiter='\t');w.writeheader();w.writerows(flat)
    (out/'measurements.json').write_text(json.dumps(measurements,indent=2),encoding='utf-8')
    render(assembly, rows, measurements, provenance, sections, context, out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--assembly',required=True);p.add_argument('--calls',required=True);p.add_argument('--context',required=True);p.add_argument('--out',required=True);a=p.parse_args();generate(a.assembly,a.calls,a.context,a.out)
