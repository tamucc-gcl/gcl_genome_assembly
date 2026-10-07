"""GitHub-readable evidence presentation; no cutting decisions are inferred."""
import json
import re


def peer_observability(t):
    if t['relationship']!='uninformative':
        return 'Both sides meet chromosome-assignment requirements'
    reasons=[]
    for side in ('left','right'):
        value=t[side]
        if value.get('qualified'):
            continue
        if not value.get('aligned_bp'):
            reason='no qualifying aligned bases after filtering; raw matches may exist'
        elif not value.get('chrom'):
            reason='aligned sequence has ambiguous or mixed chromosome assignment'
        else:
            reason='assignment to '+str(value['chrom'])+' has insufficient qualifying bases, dominance or informative bins'
        reasons.append(side+': '+reason)
    return '; '.join(reasons) or 'Chromosome relationship could not be classified'


def decision_context(m, sample):
    """Describe measured chromosome relationships; review priority never selects a cut."""
    groups = {}
    for t in m.get('chromosome_tracks', []):
        if t.get('auto_evidence') and t['sample'] != sample:
            groups.setdefault(t['sample'], set())
            if t['relationship'] in ('different_chromosomes', 'same_chromosome'):
                groups[t['sample']].add((t['relationship'], t['left']['chrom'], t['right']['chrom']))
    separate=[];same=[];unclear=[];conflict=[];pairs=set()
    for individual, values in sorted(groups.items()):
        if len(values)>1:
            conflict.append(individual)
        elif not values:
            unclear.append(individual)
        else:
            relationship,left,right=next(iter(values))
            text=individual+': '+str(left)+' → '+str(right)
            (separate if relationship=='different_chromosomes' else same).append(text)
            pairs.add(str(left)+' → '+str(right))
    for_cut='Separate chromosomes: '+'; '.join(separate) if separate else 'No informative independent chromosome evidence supporting a break'
    against='Same chromosome: '+'; '.join(same) if same else 'No opposing chromosome evidence observed'
    spans=m.get('hifi_spanning_molecules')
    if spans and m.get('hifi_informative'):
        against+='; '+str(spans)+' qualified immediate spanning molecules'
    limits=[]
    if unclear:limits.append('No informative two-sided chromosome assignment: '+', '.join(unclear))
    if conflict:limits.append('Discordant haplotype assignments within individual: '+', '.join(conflict))
    if not m.get('hifi_informative'):limits.append('Immediate HiFi assay not informative; zero spanning reads is inconclusive')
    if not m.get('hic_informative'):limits.append('Hi-C assay not calibrated/informative for a cut decision')
    if m.get('local_path_support')=='supported_grid':
        against+='; local HiFi continuity supported across the sampled transition'
    if same and not separate and not conflict:
        priority='Evidence favors retaining this sampled boundary'
    elif m.get('verified_gap') and len(separate)>=2 and not same and not conflict and not (spans and m.get('hifi_informative')) and m.get('local_path_support')!='supported_grid':
        priority='Prioritize gap-cut review'
    elif separate:
        priority='Investigate chromosome transition'+('; exact cut not localized' if not m.get('verified_gap') else '; mixed/limited evidence')
    else:
        priority='Insufficient evidence to propose a break'
    return dict(chromosome_context='; '.join(sorted(pairs)) or 'Unresolved from qualified local alignments',
        evidence_for_cut=for_cut,evidence_against_cut=against,review_priority=priority,
        evidence_limits='; '.join(limits) or 'See local measurements and controls')


def table(headers, rows):
    def cell(v):
        return str(v).replace('|', '&#124;').replace('\n', ' ')
    return ['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |'] + [
        '| ' + ' | '.join(cell(v) for v in row) + ' |' for row in rows]


def render(assembly, rows, measurements, provenance, sections, context, out):
    checksum = provenance.get('sha256') or next((r['assessment_sha256'] for r in rows), '')
    (out/'registry.json').write_text(json.dumps(dict(assembly=assembly, assessment_sha256=checksum,
        coordinate_stage='pre_finishing', candidate_count=len(rows)), indent=2), encoding='utf-8')
    md = ['# ' + assembly + ' — chimera evidence', '', '[Cohort report](../README.md) · [Cut instructions](../cut-instructions.md)', '',
          '## Assessment summary', '', str(len(rows)) + ' candidate boundaries. All selections start at NO. Evidence generation applies no cuts.', '',
          'Locations use the original pre-finishing FASTA. A proposed gap cut is a verified position for review, not approval to break.', '', '## Candidate boundaries', '']
    md += ['IDs are scoped to this assembly; use assembly plus ID when referring to a decision. Review priorities organize measured evidence and do not approve cuts.', '']
    md += table(['Candidate', 'Scaffold', 'Chromosomes left → right', 'Region to review, bp', 'Exact cut, bp', 'Review priority', 'Evidence for cutting', 'Evidence against cutting'], [
        ['['+r['id']+'](#candidate-'+r['id'].lower()+')', r['scaffold'], r['chromosome_context'], r['review_range'], r['cut_bp'] or 'Not assigned', r['review_priority'], r['evidence_for_cut'],r['evidence_against_cut']] for r in rows])
    for r, section in zip(rows, sections):
        key = r['id']; m = measurements.get(r['source_candidate_id'], {}); tracks = m.get('chromosome_tracks', [])
        md += ['', '## Candidate '+key, '', '**Scaffold:** '+r['scaffold']+'. **Region to review:** '+r['review_range']+' bp. **Exact cut:** '+str(r['cut_bp'] or 'Not assigned')+'. **Selected:** NO.', '',
               '**Why this location is a range:** '+r['localization_explanation']+'.', '',
               '**Gap interval:** '+str(r['gap_start'] or 'Unavailable')+'–'+str(r['gap_end'] or 'Unavailable')+'. **Proposed action:** '+r['action']+'.', '',
               '**Chromosomes left → right:** '+r['chromosome_context']+'.', '', '**Review priority:** '+r['review_priority']+'.', '',
               '**For cutting:** '+r['evidence_for_cut']+'.', '', '**Against cutting:** '+r['evidence_against_cut']+'.', '',
               '**Limits on the decision:** '+r['evidence_limits']+'.', '', '### Across-assembly chromosome evidence', '']
        values = []
        for t in tracks:
            role = 'Same individual' if t['sample'] == provenance.get('sample', assembly.rsplit('_hap', 1)[0]) else 'Independent comparison eligible' if t.get('auto_evidence') else 'Context only'
            values.append([t['peer'], t['sample'], role, t['left'].get('chrom') or 'Unresolved', t['right'].get('chrom') or 'Unresolved',
                '%.1f / %.1f' % (t['left']['aligned_bp']/1000, t['right']['aligned_bp']/1000),
                '%.1f%% / %.1f%%'%(100*t['left'].get('coverage',0),100*t['right'].get('coverage',0)),t['relationship'],peer_observability(t)])
        md += table(['Peer assembly', 'Individual', 'Role', 'Left chromosome', 'Right chromosome', 'Aligned kb left/right', 'Assigned coverage left/right', 'Relationship', 'Measurement adequacy / limitation'], values)
        svg = re.search(r'<svg.*?</svg>', section, re.S)
        if svg and tracks:
            filename = key+'.tracks.svg'
            (out/filename).write_text(svg.group().replace('<svg ', '<svg xmlns="http://www.w3.org/2000/svg" ', 1), encoding='utf-8')
            md += ['', '![Peer chromosome tracks]('+filename+')']
        md += ['', 'Different chromosomes means the assessed sides map to separate chromosomes in the peer, not that the peer has a fusion. Haplotypes are grouped by individual in the summary. Absence of an expected homologous match can be evidence when sequence availability and assay sensitivity are established. Failure of a qualifying alignment filter alone does not establish biological absence; the coverage and limitation columns show what was measured.', '',
               '### Local sequence and contact support', '']
        raw = m.get('hifi_raw', {})
        md += table(['Assay', 'Measurement'], [
            ['Qualified immediate HiFi spanning molecules', m.get('hifi_spanning_molecules', 'Unavailable')],
            ['Qualified HiFi flank molecules left/right', str(raw.get('left_molecules', 'Unavailable'))+' / '+str(raw.get('right_molecules', 'Unavailable'))],
            ['HiFi informative', m.get('hifi_informative', 'Unavailable')], ['Graph context', m.get('graph_status', 'Unavailable')]])
        md += ['', 'Zero spanning reads must be interpreted with flank coverage, ambiguity and interval width. Graph connectivity alone does not establish a correct join.', '']
        farther=[]
        for offset, assay in sorted(m.get('farther_hifi',{}).items(),key=lambda x:int(x[0])):
            raw=assay.get('raw',{})
            farther.append([int(offset)//1000,raw.get('left_molecules','Unavailable'),raw.get('right_molecules','Unavailable'),raw.get('spanning','Unavailable'),
                str(raw.get('left_median_depth','Unavailable'))+' / '+str(raw.get('right_median_depth','Unavailable')),assay.get('informative','Unavailable')])
        md += table(['HiFi offset kb','Left molecules','Right molecules','Spanning molecules','Median depth left/right','Flanks observable'],farther)
        md += ['', 'Observable distant flanks show reads are available on each side; no spanning reads across a long interval do not by themselves test the exact seam.', '']
        values = []
        for t in m.get('farther_contact_evidence', {}).get('trials', []):
            c = t['raw_counts']; pop = t['control_populations']
            values.append([t['library'], t['offset_bp']//1000, c.get('cross', 'Unavailable'),
                str(c.get('left_within', 'Unavailable'))+' / '+str(c.get('right_within', 'Unavailable')),
                str(pop['continuous_control'])+' / '+str(pop['gap_control']), t['informative']])
        md += table(['Library', 'Offset kb', 'Cross pairs', 'Within left/right', 'Sequence/gap controls', 'Informative'], values)
        link = '../../sequence_context/'+assembly+'.sequence_context/'
        packet = m.get('packet_interval_id')
        if packet:
            for suffix, title in [('.controls.png', 'Immediate measurements and controls'), ('.farther_contacts.png', 'Farther contact evidence')]:
                if (context/(packet+suffix)).exists():
                    md += ['', '!['+title+']('+link+packet+suffix+')']
            if (context/(packet+'.igv.xml')).exists():
                md += ['', '[IGV session: original coordinates]('+link+packet+'.igv.xml)']
        md += ['', 'Scaffolding Hi-C is corroboration, not independent validation. Small control populations and poor observability limit conclusions from weak support.', '',
               '**Decision needed:** review supporting, opposing and missing evidence before selecting an exact cut. Leave retained or unresolved rows at NO and record reviewer and rationale.', '',
               '<details><summary>Full candidate measurements</summary>', '', '```json', json.dumps(m, indent=2), '```', '', '</details>']
    md += ['', '## Source identity and measurements', '', '**Assessment SHA-256:** `'+checksum+'`.', '',
           '[Raw measurements](measurements.json) · [Source provenance](../../sequence_context/'+assembly+'.sequence_context/provenance.json)', '',
           'Track summaries use MAPQ ≥30, record identity ≥90%, ambiguity masking and coverage/dominance requirements. These exploratory measurements do not grant cutting permission.']
    (out/'report.md').write_text('\n'.join(md)+'\n', encoding='utf-8')
