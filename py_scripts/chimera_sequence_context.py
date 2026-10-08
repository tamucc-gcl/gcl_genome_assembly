#!/usr/bin/env python3
"""Published-tool diagnostic orchestration. Never emits permission to cut."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import re
from collections import defaultdict
import xml.etree.ElementTree as ET
from chimera_tracks import tracks, bins, side_summary
from chimera_blocks import chromosome_blocks, summarize_blocks, farther_hifi, farther_regions, farther_contacts
from chimera_controls import (agp, select_controls, sam_measure, peer_relationship,
                              scan_contacts, compare_controls, continuity_grid, nearby_hypotheses, continuous_controls, peer_anchor_trials)
from chimera_graph_evidence import measure_graph


def eligible(rows):
    return [r for r in rows if r.get('chromosome_member') == 'yes'
            and r.get('candidate_verdict') in ('REVIEW', 'BREAK_CANDIDATE')
            and (r.get('callable') == 'yes' or r.get('evidence_only') == 'yes')]


def spans_interval(fields, lo, hi, flank=1000):
    """A single gapped PAF record brackets an interval; NOT a proof of continuity."""
    return int(fields[2]) <= lo - flank and int(fields[3]) >= hi + flank


def run(args, output=None):
    if output:
        with open(output, 'w') as handle:
            subprocess.run([str(x) for x in args], stdout=handle, check=True)
    else:
        subprocess.run([str(x) for x in args], check=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fasta', required=True)
    p.add_argument('--calls', required=True)
    p.add_argument('--assembly', required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--peers', required=True, help='JSON list of id/sample/path records')
    p.add_argument('--reads', default='')
    p.add_argument('--bam', default='', help='Reuse a BAM bound by --bam-provenance')
    p.add_argument('--bam-provenance', default='')
    p.add_argument('--agp', default='')
    p.add_argument('--pairs', default='')
    p.add_argument('--pairs-source', default='')
    p.add_argument('--libraries', default='')
    p.add_argument('--native-graph', default='')
    p.add_argument('--motif', required=True)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    if a.reads and a.bam:
        p.error('Supply reads or an existing BAM, not both')
    if a.bam and not a.bam_provenance:
        p.error('Existing BAM reuse requires checksum provenance')
    if a.agp and not all((a.pairs, a.pairs_source, a.libraries)):
        p.error('AGP assay requires pairs, their input FASTA and library manifest')
    out = Path(a.out)
    out.mkdir()
    with open(a.calls) as handle:
        rows = eligible(list(csv.DictReader((x for x in handle if not x.startswith('#')), delimiter='\t')))
    digest = hashlib.sha256()
    with open(a.fasta, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            digest.update(chunk)
    if not rows:
        (out/'report.md').write_text('# Chimera sequence context\n\nNo eligible chromosome-scale candidate intervals.\n')
        (out/'provenance.json').write_text(json.dumps(dict(assembly=a.assembly,sample=a.sample,sha256=digest.hexdigest(),
            coordinate_stage='pre_finishing',intervals={},hifi=bool(a.reads or a.bam)),indent=2)+'\n')
        return
    if any(r.get('assembly_sha256') != digest.hexdigest() for r in rows):
        raise ValueError('Candidate assessment checksum does not match FASTA')
    # Index a local symlink; never write an index next to a published input.
    local = out/'assessment.fa'
    local.symlink_to(Path(a.fasta).resolve())
    run(['samtools', 'faidx', local])
    sizes = {f[0]: int(f[1]) for f in (x.split('\t') for x in Path(str(local)+'.fai').read_text().splitlines())}
    if sum(sizes.values()) >= 8_000_000_000:
        raise ValueError('Assessment exceeds diagnostic single-index budget')
    queries = out/'intervals.fa'
    intervals = {}
    decision_measurements = {}
    with queries.open('w') as combined:
        for i, row in enumerate(rows):
            lo, hi = int(row['transition_lo']), int(row['transition_hi'])
            chrom = row['scaffold']
            if not 0 <= lo <= hi <= sizes[chrom]:
                raise ValueError('Invalid assessment interval')
            start, end = max(0, lo-850000), min(sizes[chrom], hi+850000)
            key = 'candidate_%d' % (i+1)
            fa = out/(key+'.fa')
            run(['samtools', 'faidx', local, '%s:%d-%d' % (chrom, start+1, end)], fa)
            seq = ''.join(fa.read_text().splitlines()[1:])
            fa.write_text('>'+key+'\n'+seq+'\n')
            combined.write(fa.read_text())
            intervals[key] = dict(scaffold=chrom, start=start, end=end, lo=lo, hi=hi)
            decision_id = hashlib.sha256(json.dumps(row, sort_keys=True).encode()).hexdigest()[:20]
            intervals[key]['decision_id'] = decision_id
            decision_measurements[decision_id] = dict(assessment_sha256=digest.hexdigest(),
                                                      coordinate_stage='pre_finishing',packet_interval_id=key)
    candidate_keys = list(intervals)
    original_candidate_keys=list(candidate_keys)
    controls, links, placements = {}, {}, {}
    if a.agp:
        gaps, placements, agp_lengths = agp(a.agp)
        from chimera_controls import literal_gaps
        literal = literal_gaps(local)
        # AGP describes only this scaffolding round: embedded source gaps remain real gaps.
        for g in gaps:
            if not any(x['scaffold']==g['scaffold'] and x['lo']<=g['lo'] and x['hi']>=g['hi'] for x in literal):
                raise ValueError('AGP gap is absent from literal assessment N runs')
        gaps = literal
        with (out/'literal_gap_catalog.tsv').open('w') as catalog:
            writer=csv.DictWriter(catalog,fieldnames=['scaffold','lo','hi'],delimiter='\t')
            writer.writeheader();writer.writerows(gaps)
        if agp_lengths != sizes:
            raise ValueError('AGP object dictionary does not match assessment FASTA')
        hypotheses=nearby_hypotheses(intervals,gaps)
        for key,interval in hypotheses.items():
            decision_id=hashlib.sha256(json.dumps([a.assembly,digest.hexdigest(),interval['scaffold'],interval['lo'],interval['hi']]).encode()).hexdigest()[:20]
            interval['decision_id']=decision_id
            decision_measurements[decision_id]=dict(assembly=a.assembly,scaffold=interval['scaffold'],assessment_sha256=digest.hexdigest(),
                coordinate_stage='pre_finishing',review_only=True,source_candidates=interval['source_candidates'],packet_interval_id=key,
                gap_start=interval['lo'],gap_end=interval['hi'],cut_bp=interval['hi'])
        intervals.update(hypotheses);candidate_keys+=list(hypotheses)
        controls, links = select_controls(intervals, gaps, sizes, count=48)
        continuous,continuous_links=continuous_controls(intervals,gaps,sizes)
        controls.update(continuous)
        for key in intervals:links[key]=links.get(key,[])+continuous_links.get(key,[])
        for key in candidate_keys:
            interval=intervals[key]
            interval['transition_lo'],interval['transition_hi']=interval['lo'],interval['hi']
            if 'gap' in interval:
                # Absence/contact assays test the literal gap, never a broad core.
                interval['lo'],interval['hi']=interval['gap']['lo'],interval['gap']['hi']
        with queries.open('a') as combined:
            for key, interval in dict(controls,**hypotheses).items():
                start, end = max(0, interval['lo']-850000), min(sizes[interval['scaffold']], interval['hi']+850000)
                interval.update(start=start, end=end)
                fa = out/(key+'.fa')
                run(['samtools', 'faidx', local, '%s:%d-%d' % (interval['scaffold'], start+1, end)], fa)
                seq = ''.join(fa.read_text().splitlines()[1:])
                fa.write_text('>'+key+'\n'+seq+'\n'); combined.write(fa.read_text())
            intervals.update(controls)
        for key, interval in intervals.items():
            if 'gap' in interval:
                seq=''.join((out/(key+'.fa')).read_text().splitlines()[1:])
                g=interval['gap']
                if set(seq[g['lo']-interval['start']:g['hi']-interval['start']].upper()) != {'N'}:
                    raise ValueError('AGP gap is not literal Ns in assessment')
        for key in candidate_keys:
            interval=intervals[key];measurement=decision_measurements[interval['decision_id']]
            measurement.update(verified_gap='gap' in interval,assessment_scaffold_length=sizes[interval['scaffold']])
            if 'gap' in interval:measurement.update(gap_start=interval['gap']['lo'],gap_end=interval['gap']['hi'],cut_bp=interval['gap']['hi'])
        # Also expose nearby gaps, without silently snapping the proposed boundary.
        nearby = [dict(candidate=key, **g, relationship='within_transition' if interval['transition_lo']<=g['lo']<g['hi']<=interval['transition_hi'] else 'nearby_not_localized')
                  for key in original_candidate_keys for interval in [intervals[key]] for g in gaps
                  if g['scaffold']==interval['scaffold'] and g['lo']<interval['transition_hi']+250000 and g['hi']>interval['transition_lo']-250000]
        if nearby:
            with (out/'gap_hypotheses.tsv').open('w') as handle:
                w=csv.DictWriter(handle,fieldnames=list(nearby[0]),delimiter='\t');w.writeheader();w.writerows(nearby)
    run(['tidk', 'search', '--string', a.motif, '--window', '100', '--output', 'intervals', '--dir', out/'tidk', queries])
    peers = json.loads(a.peers)
    summaries = {key: dict(same=set(), other=set()) for key in intervals}
    peer_assays = defaultdict(lambda: defaultdict(set))
    peer_rows=[]
    anchor_rows=[]
    block_rows=defaultdict(list)
    track_rows=defaultdict(list)
    with (out/'peer_context.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['candidate', 'peer', 'sample', 'relationship', 'single_record_brackets_interval', 'mapq', 'target', 'target_start', 'target_end', 'target_length'])
        for n, peer in enumerate(peers):
            peer_local = out/('peer_%d.fa' % n)
            peer_local.symlink_to(Path(peer['path']).resolve())
            run(['samtools', 'faidx', peer_local])
            length = sum(int(line.split('\t')[1]) for line in Path(str(peer_local)+'.fai').read_text().splitlines())
            if length >= 8_000_000_000:
                raise ValueError('Peer exceeds diagnostic single-index budget')
            paf = out/('peer_%d.paf' % n)
            run(['minimap2', '-x', 'asm5', '-I', '8G', '-c', '--secondary=yes', '-N', '50', '-p', '0.5', '-t', a.threads, peer['path'], queries], paf)
            # Report every alignment, including ambiguity; a missing row is not negative support.
            grouped=defaultdict(list)
            for line in paf.read_text().splitlines():
                f = line.split('\t')
                grouped[f[0]].append(f)
                interval = intervals[f[0]]
                bracket = spans_interval(f, interval['lo']-interval['start'], interval['hi']-interval['start'])
                if bracket and int(f[11]) >= 20:
                    relationship = 'same' if peer['sample'] == a.sample else 'other'
                    summaries[f[0]][relationship].add(peer['sample'])
                writer.writerow([f[0], peer['id'], peer['sample'], 'same_individual' if peer['sample'] == a.sample else 'other_individual', 'yes' if bracket else 'no', f[11], f[5], f[7], f[8], f[6]])
            for key, interval in intervals.items():
                (relationship,left,right),trials=peer_anchor_trials(grouped[key],interval['lo']-interval['start'],interval['hi']-interval['start'],interval['end']-interval['start'])
                anchor_rows.extend(dict(candidate=key,peer=peer['id'],**trial) for trial in trials)
                block_rows[key].append(chromosome_blocks(grouped[key],interval,peer,gaps if a.agp else []))
                if 'decision_id' in interval:
                    track=tracks(grouped[key],peer.get('chromosome_labels',{}))
                    tiles=bins(track,interval['end']-interval['start'])
                    track_left=side_summary(track,tiles,0,interval['lo']-interval['start'])
                    track_right=side_summary(track,tiles,interval['hi']-interval['start'],interval['end']-interval['start'])
                    relation=('different_chromosomes' if track_left['chrom']!=track_right['chrom'] else 'same_chromosome') if track_left['qualified'] and track_right['qualified'] else 'uninformative'
                    track_rows[key].append(dict(peer=peer['id'],sample=peer['sample'],auto_evidence=peer.get('auto_evidence') is True,
                        relationship=relation,left=track_left,right=track_right,bins=tiles,
                        bridge_segments=[b for b in track if b['lo']<interval['hi']-interval['start']+250000 and b['hi']>interval['lo']-interval['start']-250000]))
                labels=peer.get('chromosome_labels',{})
                left_chrom=labels.get(left[5]) if left else None
                right_chrom=labels.get(right[5]) if right else None
                if relationship=='separate_scaffolds':
                    relationship=('different_chromosomes' if left_chrom and right_chrom and left_chrom!=right_chrom else
                                  'fragmented_same_chromosome' if left_chrom and left_chrom==right_chrom else 'uninformative_chromosome_identity')
                if relationship=='continuous_context' and not left_chrom:
                    relationship='uninformative_chromosome_identity'
                if peer['sample'] != a.sample and peer.get('auto_evidence') is True and relationship != 'uninformative':
                    peer_assays[key][peer['sample']].add(relationship)
                peer_rows.append(dict(candidate=key,peer=peer['id'],sample=peer['sample'],relationship=relationship,auto_evidence=peer.get('auto_evidence') is True,
                    left_chrom=left_chrom or '.',right_chrom=right_chrom or '.',
                    left_target=left[5] if left else '.',right_target=right[5] if right else '.',
                    left_target_start=left[7] if left else '.',left_target_end=left[8] if left else '.',
                    right_target_start=right[7] if right else '.',right_target_end=right[8] if right else '.',
                    left_strand=left[4] if left else '.',right_strand=right[4] if right else '.'))
            peer_local.unlink()
    if anchor_rows:
        with (out/'peer_anchor_trials.tsv').open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(anchor_rows[0]),delimiter='\t');writer.writeheader();writer.writerows(anchor_rows)
    flat_blocks=[]
    for key,assays in block_rows.items():
        for assay in assays:
            for trial in assay['trials']:
                flat_blocks.append(dict(candidate=key,peer=assay['peer'],sample=assay['sample'],auto_evidence=assay['auto_evidence'],
                    qualified=assay['qualified'],unique_gap_localization=assay['unique_gap_localization'],**trial))
    if flat_blocks:
        with (out/'chromosome_blocks.tsv').open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(flat_blocks[0]),delimiter='\t');writer.writeheader();writer.writerows(flat_blocks)
    if peer_rows:
        with (out/'peer_boundary_assays.tsv').open('w') as handle:
            w=csv.DictWriter(handle,fieldnames=list(peer_rows[0]),delimiter='\t');w.writeheader();w.writerows(peer_rows)
    if a.reads or a.bam:
        if a.bam:
            provenance=json.loads(Path(a.bam_provenance).read_text())
            if provenance.get('sha256', provenance.get('assessment_sha256')) != digest.hexdigest() or provenance.get('coordinate_stage') != 'pre_finishing':
                raise ValueError('BAM provenance does not match assessment FASTA')
            bam=Path(a.bam).resolve()
            dictionary=out/'bam_dictionary.tsv'
            run(['samtools','idxstats',bam],dictionary)
            bam_sizes={f[0]:int(f[1]) for f in (line.split('\t') for line in dictionary.read_text().splitlines()) if f[0]!='*'}
            if bam_sizes != sizes: raise ValueError('BAM dictionary differs from assessment')
        else:
            sam = out/'mapping.sam'
            run(['minimap2', '-ax', 'map-hifi', '-I', '8G', '--secondary=yes', '-N', '50', '-p', '0.5', '-t', a.threads, local, a.reads], sam)
            bam = out/'hifi.bam'
            run(['samtools', 'sort', '-@', '2', '-m', '2G', '-o', bam, sam])
            sam.unlink()
            run(['samtools', 'index', bam])
        run(['samtools', 'quickcheck', bam])
        for key, interval in intervals.items():
            region = '%s:%d-%d' % (interval['scaffold'], interval['start']+1, interval['end'])
            run(['samtools', 'view', '-h', bam, region], out/(key+'.sam'))
            interval['hifi']=sam_measure(out/(key+'.sam'),interval['lo'],interval['hi'])
            interval['farther_hifi']=farther_hifi(out/(key+'.sam'),interval)
            if 'decision_id' in interval:
                decision_measurements[interval['decision_id']]['hifi_spanning_molecules'] = interval['hifi']['spanning']
                if 'gap' not in interval:
                    grid=continuity_grid(out/(key+'.sam'),interval['lo'],interval['hi'])
                    decision_measurements[interval['decision_id']].update(local_path_support='supported_grid' if grid['minimum_molecules']>=2 else 'unresolved',continuity_grid=grid)
            region='%s:%d-%d'%(interval['scaffold'],max(1,interval['lo']-100000+1),min(sizes[interval['scaffold']],interval['hi']+100000))
            run(['samtools', 'depth', '-aa', '-q', '0', '-Q', '20', '-G', '0xF04', '-r', region, bam], out/(key+'.depth.tsv'))
    for key, interval in intervals.items():
        interval['key']=key
        interval['peer_continuous_individuals']=sum(value=={'continuous_context'} for value in peer_assays[key].values())
        interval['peer_separate_individuals']=sum(value=={'different_chromosomes'} for value in peer_assays[key].values())
    for key in candidate_keys:
        decision_measurements[intervals[key]['decision_id']]['chromosome_tracks']=track_rows[key]
    graph_assays=measure_graph(a.native_graph,out,queries,intervals,run,a.threads) if a.native_graph else {}
    for key in candidate_keys:
        decision_measurements[intervals[key]['decision_id']].update(graph_assays.get(key,dict(graph_status='native_graph_unavailable',graph_contradiction=None)))
    if a.agp and (a.reads or a.bam):
        source=out/'pairs_input.fa'
        source.symlink_to(Path(a.pairs_source).resolve());run(['samtools','faidx',source])
        source_sizes={f[0]:int(f[1]) for f in (line.split('\t') for line in Path(str(source)+'.fai').read_text().splitlines())}
        pair_sizes={}
        opener=__import__('gzip').open if a.pairs.endswith('.gz') else open
        with opener(a.pairs,'rt') as handle:
            for line in handle:
                if not line.startswith('#'): break
                if line.startswith('#chromsize:'):
                    f=line.split();pair_sizes[f[1]]=int(f[2])
        if pair_sizes != source_sizes: raise ValueError('Hi-C pairs dictionary differs from exact input FASTA')
        with open(a.libraries) as handle: libraries=list(csv.DictReader(handle,delimiter='\t'))
        far_regions=farther_regions(intervals,sizes)
        counts,totals,audit=scan_contacts(a.pairs,placements,dict(intervals,**far_regions),libraries)
        raw=[]
        for (library,key),value in counts.items():
            raw.append(dict(library=library,candidate=key,**{field:value[field] for field in ('cross','left_within','right_within','left_ends','right_ends')}))
        if raw:
            with (out/'contact_controls.tsv').open('w') as handle:
                w=csv.DictWriter(handle,fieldnames=list(raw[0]),delimiter='\t');w.writeheader();w.writerows(raw)
        for key in candidate_keys:
            interval=intervals[key];measurement=decision_measurements[interval['decision_id']]
            comparison=compare_controls(interval,[intervals[c] for c in links.get(key,[])],counts,{r['library_id'] for r in libraries})
            block=summarize_blocks(block_rows[key],a.sample)
            far=farther_contacts(interval,[intervals[c] for c in links.get(key,[])],counts,{r['library_id'] for r in libraries},far_regions)
            own_pairs={tuple(r['chromosome_pair']) for r in block_rows[key] if r['sample']==a.sample and r['qualified']}
            measurement.update(chromosome_blocks=block,farther_contact_evidence=far,farther_hifi=interval['farther_hifi'],
                haplotype_block_conflict=bool(own_pairs and block['chromosome_pair'] and own_pairs!={tuple(block['chromosome_pair'])}),
                repeat_obscured_localization=block['localized'] and 'gap' in interval)
            measurement.update(comparison, independent_discordant_individuals=interval['peer_separate_individuals'],
                verified_gap='gap' in interval, alternative_placements_checked=True,
                control_ids=links.get(key,[]), hifi_raw=interval['hifi'],
                alternative_placement_scope='Assessment BAM MAPQ and competing peer anchors; bounded emitted alternatives, no proof of haplotype-specific uniqueness')
            measurement['assessment_scaffold_length']=sizes[interval['scaffold']]
            if 'gap' in interval:
                measurement.update(gap_start=interval['gap']['lo'],gap_end=interval['gap']['hi'],cut_bp=interval['gap']['hi'])
        source.unlink()
        (out/'contact_audit.json').write_text(json.dumps(dict(audit=dict(audit),library_totals=dict(totals),
            interpretation='Scaffolding reads are corroboration, not independent validation; thresholds are a conservative heuristic.'),indent=2)+'\n')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for key in candidate_keys:
            interval=intervals[key]
            ids=[key]+links.get(key,[])
            library_ids=sorted({r['library_id'] for r in libraries})
            fig,axes=plt.subplots(len(library_ids)+1,1,figsize=(max(12,len(ids)*.35),3*(len(library_ids)+1)),squeeze=False)
            for index,library in enumerate(library_ids):
                ax=axes[index,0]
                for offset,field,color in ((-.25,'cross','tab:red'),(0,'left_within','tab:blue'),(.25,'right_within','tab:green')):
                    ax.bar([n+offset for n in range(len(ids))],[counts[library,c][field] for c in ids],width=.25,label=field,color=color)
                ax.set_xticks(range(len(ids)),['candidate']+[('gap ' if intervals[c]['role']=='gap_control' else 'sequence ')+str(n+1) for n,c in enumerate(ids[1:])],rotation=30)
                ax.set_ylabel('UU pairs in 250 kb flanks');ax.set_title(library);ax.legend()
            ax=axes[-1,0]
            for offset,field in ((-.25,'spanning'),(0,'left_molecules'),(.25,'right_molecules')):
                ax.bar([n+offset for n in range(len(ids))],[intervals[c]['hifi'][field] for c in ids],width=.25,label=field)
            ax.set_xticks(range(len(ids)),['candidate']+[('gap ' if intervals[c]['role']=='gap_control' else 'sequence ')+str(n+1) for n,c in enumerate(ids[1:])],rotation=30)
            ax.set_ylabel('Distinct MAPQ30 molecules');ax.legend()
            fig.suptitle('%s %s:%d-%d — raw candidate/control assays, not fusion probabilities' % (a.assembly,interval['scaffold'],interval['lo'],interval['hi']))
            fig.tight_layout();fig.savefig(out/(key+'.controls.png'),dpi=150);plt.close(fig)
        for key in candidate_keys:
            measurement=decision_measurements[intervals[key]['decision_id']]
            trials=measurement.get('farther_contact_evidence',{}).get('trials',[])
            fig,axes=plt.subplots(len(library_ids),1,figsize=(10,3*len(library_ids)),squeeze=False)
            for index,library in enumerate(library_ids):
                values=[t for t in trials if t['library']==library]
                ax=axes[index,0]
                ax.bar([n-.2 for n in range(len(values))],[t['upper_count_allowance_ratio'] for t in values],width=.4,label='candidate ratio + count allowance')
                ax.bar([n+.2 for n in range(len(values))],[.1*t['minimum_control_ratio'] for t in values],width=.4,label='10% of minimum matched-control ratio')
                ax.set_xticks(range(len(values)),[str(t['offset_bp']//1000)+' kb; '+('informative' if t['informative'] else 'uninformative') for t in values],rotation=15)
                ax.set_ylabel('Normalized cross-flank contact');ax.set_title(library);ax.legend()
            fig.suptitle('%s %s - farther-flank assays; uninformative points cannot authorize a cut'%(a.assembly,key))
            fig.tight_layout();fig.savefig(out/(key+'.farther_contacts.png'),dpi=150);plt.close(fig)
    if a.reads or a.bam:
        for key in candidate_keys:
            interval=intervals[key]
            session=ET.Element('Session',genome=str(Path(a.fasta).resolve()),locus='%s:%d-%d' % (interval['scaffold'],interval['start']+1,interval['end']),version='8')
            resources=ET.SubElement(session,'Resources');ET.SubElement(resources,'Resource',path=str(bam))
            ET.ElementTree(session).write(out/(key+'.igv.xml'),encoding='utf-8',xml_declaration=True)
    (out/'control_registry.json').write_text(json.dumps(dict(intervals=intervals,control_links=links),indent=2)+'\n')
    (out/'decision_measurements.json').write_text(json.dumps(decision_measurements, indent=2)+'\n')
    for tool in ('samtools', 'minimap2', 'tidk'):
        run([tool, '--version'], out/(tool+'.version.txt'))
    local.unlink()
    (out/'provenance.json').write_text(json.dumps(dict(assembly=a.assembly, sample=a.sample, sha256=digest.hexdigest(), coordinate_stage='pre_finishing', intervals=intervals, peers=peers, hifi=bool(a.reads or a.bam), reused_bam=a.bam or None), indent=2)+'\n')
    (out/'report.md').write_text(
        '# Chimera sequence context: '+a.assembly+'\n\n'
        +str(len(candidate_keys))+' candidate intervals and '+str(len(controls))+' comparison sites; '+str(len(peers))+' same-species peer assemblies.\n\n'
        'HiFi mapping: '+('available' if a.reads or a.bam else 'not requested or no HiFi reads')+'.\n\n'
        'See peer_context.tsv for direct interval alignments and provenance.json for coordinate offsets. '
        'A single gapped alignment bracketing an interval is contextual support, not proof of basewise continuity. '
        'Same-individual haplotypes are not independent votes. No alignment is not evidence of a technical error.\n\n'
        'tidk windows describe motif counts, not repeat-tract lengths or chromosome fusions. '
        'HiFi evidence requires inspection of anchors, repeats and alternative placements. '
        'No diagnostic here authorizes a cut; candidate and cutting decisions remain separate.\n')
    with (out/'report.md').open('a') as report:
        report.write('\n## Boundary/control assays\n\nSee decision_measurements.json and control_registry.json for complete measured values and literal coordinates. '
            'Automatic gap eligibility requires at least three usable independent individuals and five matched controls from each separately qualified population; both libraries must be informative and show support loss. '
            'These thresholds are heuristic and require validation. Native graph placement status and its limits are recorded in the packet. '
            'IGV sessions point to the exact pre-finishing reference/BAM; paths need updating if the files are copied to a workstation.\n\n')
        for key in candidate_keys:
            report.write('### '+key+'\n\n')
            if (out/(key+'.controls.png')).exists():report.write('!['+key+' candidate/control measurements]('+key+'.controls.png)\n\n')
            if (out/(key+'.farther_contacts.png')).exists():report.write('!['+key+' farther flank contacts]('+key+'.farther_contacts.png)\n\n')
            if (out/(key+'.igv.xml')).exists():report.write('[IGV session]('+key+'.igv.xml)\n\n')
        report.write('\n<details>\n<summary>Interval-level peer context</summary>\n\n'
                     'Counts below require MAPQ >=20 and a single alignment with 1-kb flanks. '
                     'They count distinct individuals, not alignment records, and are not cutting votes.\n\n'
                     '| Candidate | Scaffold interval (0-based, half-open) | Same individual observed | Other individuals observed |\n'
                     '|---|---|---|---:|\n')
        for key, interval in intervals.items():
            summary = summaries[key]
            report.write('| %s | %s:%d-%d | %s | %d |\n' %
                         (key, interval['scaffold'], interval['lo'], interval['hi'],
                          'yes' if summary['same'] else 'not observed', len(summary['other'])))
        report.write('\nZero or not observed means no qualifying record, not demonstrated discontinuity.\n\n</details>\n')


if __name__ == '__main__':
    main()
