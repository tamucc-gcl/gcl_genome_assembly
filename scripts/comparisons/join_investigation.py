#!/usr/bin/env python3
"""Read-only junction investigation; inputs and hypotheses are supplied in a manifest.

This deliberately does not change production discovery, scoring, or assemblies.
All intervals are zero-based, half-open; a probe is the number of bases before a cut.
"""
import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import statistics
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'py_scripts'))
from chimera_controls import agp, literal_gaps, molecule_chains, qualified_molecule
from chimera_tracks import blocks, LabelIndex, position
from misassembly_core import read_paf, independent_labels, sha

from misassembly_junction import (table, write, save, run, windows, local_hits, union, overlap, consensus_masks, local_tracks, labels, sam_records, span_record, tag, alignment_summary, arm_coverage, minimum_anchor_span, hifi, library_id, window_counts, aggregate_matrix, mask_bins, hic, COMMANDS, credible_bridge)


def command_args(task, executable):
    """Read the actual cached command, never infer an input from a similar filename."""
    tokens = shlex.split((task / '.command.sh').read_text())
    index = next((i for i, t in enumerate(tokens) if Path(t).name == executable), None)
    if index is None:
        raise ValueError(f'{executable} absent from {task}/.command.sh')
    result = {}
    for i in range(index + 1, len(tokens)-1):
        if tokens[i].startswith('--'):
            result[tokens[i][2:]] = tokens[i+1]
    return result


def required(path):
    path = Path(path).resolve()
    if not path.is_file():
        raise ValueError(f'Required cached input missing: {path}')
    return str(path)


def stage_commands(task):
    """Read only nxf_stage, excluding similarly named output-copy commands."""
    wrapper=task/'.command.run'
    if not wrapper.is_file():return []
    commands=[];inside=False
    for line in wrapper.read_text().splitlines():
        if re.match(r'^\s*nxf_stage\(\)\s*\{\s*$',line):inside=True;continue
        if inside and re.match(r'^\s*}\s*$',line):break
        if inside and line.strip() and not line.lstrip().startswith('#'):commands.append(line.strip())
    return commands


def shell_commands(line):
    """Tokenize literal simple commands without executing mkdir, copies or links.

    Nextflow prefixes nested stage names with `mkdir -p DIR && ...` on one
    line. Shell punctuation must be tokenized before selecting cp/ln commands;
    punctuation inside a quoted filename remains part of that filename.
    """
    lexer=shlex.shlex(line,posix=True,punctuation_chars=';&|')
    lexer.whitespace_split=True;lexer.commenters='#'
    command=[]
    for token in lexer:
        if token in ('&&','||',';','&','|'):
            if command:yield command;command=[]
        else:command.append(token)
    if command:yield command


def staged_links(task):
    """Read literal input source/destination mappings for copy and link modes."""
    result=[]
    for line in stage_commands(task):
        try:commands=list(shell_commands(line))
        except ValueError:continue
        for tokens in commands:
            if tokens[0] not in ('cp','ln'):continue
            operands=[t for t in tokens[1:] if not t.startswith('-')]
            if len(operands)!=2:continue
            source,destination=operands
            # Nextflow writes concrete paths. Never expand shell expressions.
            if any(c in source+destination for c in ('$','`','\n')):continue
            result.append((Path(os.path.abspath(task/destination)),Path(os.path.abspath(task/source))))
    return result


class RetainedInputError(ValueError):
    def __init__(self,requested,attempts,wrappers):
        self.details=dict(requested=str(requested),attempted_paths=attempts,wrappers=wrappers)
        super().__init__(f'Required retained input unavailable: {requested}. '
                         'No existing source was recoverable from recorded Nextflow input staging.')


def retained_input(task,value,hint=None):
    """Recover a cleaned scratch-stage input through its exact recorded link chain.

    A staged directory link can supply a file beneath it. Upstream staging links
    can also have been removed; follow their wrappers recursively. No work-tree
    search by basename, copied input, or recreated symlink is used.
    """
    requested=Path(os.path.abspath(task/value));seen=set();attempts=[];wrappers={}
    def resolve(path):
        path=Path(os.path.abspath(path))
        if path.is_file():return str(path.resolve())
        if path in seen or len(seen)>=64:return None
        seen.add(path);attempts.append(str(path))
        for parent in path.parents:
            matches=[]
            if (parent/'.command.run').is_file() and str(parent) not in wrappers:
                wrappers[str(parent)]=dict(path=str(parent/'.command.run'),commands=stage_commands(parent),
                    mappings=[dict(destination=str(d),source=str(s),source_exists=s.is_file() or s.is_dir())
                              for d,s in staged_links(parent)])
            for destination,source in staged_links(parent):
                try:rest=path.relative_to(destination)
                except ValueError:continue
                matches.append((len(destination.parts),source/rest))
            for _,source in sorted(matches,key=lambda item:item[0],reverse=True):
                found=resolve(source)
                if found:return found
        return None
    # Absolute retained paths already recorded in the evidence take precedence
    # over transient staging names; downstream checksum checks still apply.
    if hint:
        found=resolve(Path(hint))
        if found:return found
    found=resolve(requested)
    if found:return found
    raise RetainedInputError(requested,attempts,list(wrappers.values()))


def cached_task(evidence):
    # mapping_manifest provenance is a cached <assembly>.evidence/evidence.json.
    return Path(evidence).parent.parent


def prepare(a):
    """Resolve exact retained task inputs from the published packet and its cache."""
    a.out.mkdir(parents=True, exist_ok=True)
    source = a.source.resolve()
    mappings = json.loads((source/'mapping_manifest.json').read_text())
    selected = table(a.focus)
    if len({(r['assembly'], r['event']) for r in selected}) != len(selected):
        raise ValueError('Duplicate investigation target')
    grouped = defaultdict(list)
    for row in selected: grouped[row['assembly']].append(row)
    jobs = [];failures=[];resolution=dict(assemblies=[])
    def resolve_input(audit,key,task,value,hint=None):
        try:
            path=retained_input(task,value,hint)
            audit['inputs'][key]=dict(status='resolved',path=path)
            return path
        except RetainedInputError as error:
            audit['inputs'][key]=dict(status='unavailable',**error.details)
            failures.append(audit['assembly']+': '+key)
            return None
        finally:save(a.out/'input_resolution.json',resolution)
    for assembly, rows in grouped.items():
        published = source/'assemblies'/assembly/'evidence.json'
        data = json.loads(published.read_text())
        task = cached_task(mappings[assembly]['provenance'])
        args = command_args(task, 'misassembly_assess.py')
        audit=dict(assembly=assembly,assessment_task=str(task),inputs={});resolution['assemblies'].append(audit)
        hints={'fasta':data.get('assessment_fasta') or mappings[assembly].get('fasta')}
        paths = {key: resolve_input(audit,key,task,value,hints.get(key)) for key, value in args.items()
                 if key in ('fasta', 'discovery', 'mapping', 'agp', 'pairs', 'source', 'libraries')}
        bam=resolve_input(audit,'bam',task,mappings[assembly]['bam'])
        index=resolve_input(audit,'index',task,mappings[assembly]['index'])
        for key in ('fasta', 'discovery', 'mapping', 'agp', 'pairs', 'source', 'libraries'):
            if key not in paths: raise ValueError(f'{assembly}: {key} unavailable in retained task')
        for key, expected in (('fasta', data['assessment_sha256']),
                              ('agp', data['hic_audit']['agp_sha256']),
                              ('source', data['hic_audit']['source_sha256'])):
            if paths[key] is not None and sha(paths[key]) != expected: raise ValueError(f'{assembly}: changed {key}')
        if data.get('coordinate_stage') != 'pre_finishing':
            raise ValueError('Expected explicit pre_finishing evidence coordinates')
        if paths['discovery'] is None:continue
        # The staged discovery path resolves back to its actual Nextflow task.
        discovery_task = Path(paths['discovery']).parent.parent
        discovery_args = command_args(discovery_task, 'misassembly_discover.py')
        catalog = resolve_input(audit,'catalog',discovery_task,discovery_args['catalog'])
        alignments = json.loads(discovery_args['alignments'])
        for pair in alignments:
            pair['path'] = resolve_input(audit,'alignment:'+pair['a']+':'+pair['b'],discovery_task,pair['path'])
        if (catalog is None or bam is None or index is None or
                any(path is None for path in paths.values()) or any(pair['path'] is None for pair in alignments)):continue
        if { (r['a'], r['b'], r['sha256']) for r in data['alignments'] } != {
                (r['a'], r['b'], sha(r['path'])) for r in alignments }:
            raise ValueError('Cohort alignment provenance differs from reviewed evidence')
        mapping = json.loads(Path(paths['mapping']).read_text())
        if (mapping['assessment_sha256'] != data['assessment_sha256'] or
                mapping.get('assembly')!=assembly or mapping.get('sample')!=data['sample'] or
                mapping.get('coordinate_stage')!='pre_finishing'):
            raise ValueError('BAM provenance differs from assessment sequence')
        if 'mapping' in data and Path(mappings[assembly]['bam']).resolve()!=Path(data['mapping']['bam']).resolve():
            raise ValueError('Published BAM manifest differs from evidence provenance')
        run(['samtools', 'quickcheck', bam])
        sizes = {f[0]:int(f[1]) for f in (line.split('\t') for line in
                 subprocess.check_output(['samtools','idxstats',bam],text=True).splitlines()) if f[0]!='*'}
        fasta_sizes = {}
        # Do not write indexes beside immutable cache inputs.
        with open(paths['fasta']) as handle:
            name = None
            for line in handle:
                if line.startswith('>'): name=line[1:].split()[0]; fasta_sizes[name]=0
                else: fasta_sizes[name] += len(line.strip())
        if sizes != fasta_sizes: raise ValueError('BAM/FASTA sequence dictionaries differ')
        events=[]
        for row in rows:
            event=next(e for e in data['events'] if e['id']==row['event'])
            sites,lo,hi=windows(event,sizes[event['scaffold']])
            events.append(dict(event=event['id'],role=row['role'],scaffold=event['scaffold'],
                               left=event['left'],right=event['right'],lo=lo,hi=hi,sites=sites,
                               transition_lo=event['lo'],transition_hi=event['hi']))
        final_map = required(a.final_dir/(assembly+'.name_map.tsv'))
        # Names after finishing are context, not an exact coordinate lift: end trimming
        # and extension can change positions. Never silently add a parent offset.
        jobs.append(dict(assembly=assembly,sample=data['sample'],sizes=sizes,events=events,
                         evidence=required(published),catalog=catalog,alignments=alignments,
                         bam=bam,index=index,final_map=final_map,paths=paths,
                         source_evidence_sha256=sha(published),catalog_sha256=sha(catalog),final_map_sha256=sha(final_map),
                         input_stats={str(p):dict(size=Path(p).stat().st_size,mtime_ns=Path(p).stat().st_mtime_ns)
                                      for p in [bam,index,*paths.values()]},
                         assessment_sha256=data['assessment_sha256'],
                         coordinate_stage='pre_finishing',
                         hic_expected=data['hic_audit']['audit']['input_pairs']))
    if failures:
        raise ValueError('Unresolved investigation inputs: '+', '.join(failures)+'. '
                         f'Details: {a.out / "input_resolution.json"} (included in the review archive).')
    if len({j['sample'] for j in jobs}) != 1:
        raise ValueError('This competitive read investigation expects haplotypes of one individual')
    for job in jobs:
        job['competitive_fastas']={j['assembly']:j['paths']['fasta'] for j in jobs}
        save(a.out/(job['assembly']+'.json'),job)
    save(a.out/'manifest.json',dict(jobs=jobs,source=str(source),focus=selected,
                                  coordinate_convention='zero-based half-open; cut_bp bases before cut'))


def shade_masks(ax,mask,colors,alpha):
    """One collection per chromosome, never one plot object per alignment block."""
    from matplotlib.collections import PolyCollection
    grouped=defaultdict(list)
    for b in mask:
        lo=b['lo']/1e6;hi=b['hi']/1e6
        grouped[colors.get(b['chrom'],'grey')].append([(lo,0),(hi,0),(hi,1),(lo,1)])
    for color,vertices in grouped.items():
        ax.add_collection(PolyCollection(vertices,facecolors=color,edgecolors='none',alpha=alpha,
                                        transform=ax.get_xaxis_transform()),autolim=False)


def report(a):
    """One review packet, with descriptive evidence and explicit unselected actions."""
    import shutil
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    manifest=json.loads(a.manifest.read_text());a.out.mkdir(exist_ok=False)
    lines=['# Four-join investigation','',
           'Read-only evidence review. **No cuts are selected or applied.**', '',
           'Coordinates below refer to the checksum-verified pre-finishing assessment FASTA. '
           'Final chromosome names are context only: finishing can shift coordinates, so they are not an exact coordinate lift.', '',
           'Blue and orange identify the two chromosome arms within each case. '
           'HiFi spanning molecules are long reads; Hi-C contacts are pairs of short reads.', '',
           '## Cases','', '| Assembly / case | Chromosomes | Role | Evidence |', '|---|---|---|---|']
    case_lines=[];actions=[];summaries=[]
    for job in manifest['jobs']:
        assembly=job['assembly'];dest=a.out/assembly;dest.mkdir()
        for branch in ('labels','hifi','hic'):
            folder=next(p for p in a.packets if p.name==assembly+'.'+branch)
            shutil.copytree(folder,dest/branch)
        hifi=table(dest/'hifi'/'spanning_profile.tsv');original=table(dest/'hifi'/'original_bridges.tsv')
        bridges=table(dest/'hifi'/'competitive_bridges.tsv');contacts=table(dest/'hic'/'contact_profile.tsv')
        reach=table(dest/'hifi'/'unique_anchor_reach.tsv')
        tracks_data=table(dest/'labels'/'tracks.tsv');anchors=table(dest/'hic'/'unique_arm_contacts.tsv')
        masks=json.loads((dest/'labels'/'masks.json').read_text());maps=json.loads((dest/'hic'/'summary.json').read_text())['maps']
        names=table(job['final_map']);write(dest/'current_names.tsv',names,list(names[0]))
        for e in job['events']:
            key=assembly+'-'+e['event'];anchor=key.lower().replace('_','-')
            role='Target' if e['role']=='target' else 'Previously reviewed break (comparison)'
            lines.append(f"| {assembly} / {e['event']} | {e['left']} â†’ {e['right']} | {role} | [Review](#{anchor}) |")
            case_lines += ['',f'<a id="{anchor}"></a>',f"## {assembly} / {e['event']}: {e['left']} â†’ {e['right']}",'',
                           '### What this batch can establish','',
                           'Compare the chromosome transition with every measured support trough. '
                           'A small alignment island beyond a trough does not veto a cut; inspect its unique-locus support. '
                           'A trough away from the arm transition may instead be a separate local assembly problem.', '',
                           f"Assessment scaffold: `{e['scaffold']}`; transition interval **{e['transition_lo']:,}â€“{e['transition_hi']:,} bp**.", '']
            matches=[r['new_name'] for r in names if r['old_name']==e['scaffold'] or r['old_name'].startswith(e['scaffold']+'_sub_')]
            case_lines += ['Current named descendants: '+(', '.join('`'+v+'`' for v in matches) or 'No matching descendant in supplied name map.')+'.','',
                           '### Sites and molecule evidence','',
                           '| Site | Position, bp | Original HiFi â‰¥20 / â‰¥30 | HiFi spans after checking alternatives | Spans joining â‰¥1 kb unique sequence from both chromosomes |',
                           '|---|---:|---:|---:|---:|']
            for s in e['sites']:
                old=next(r for r in original if r['site']==s['id'] and r['anchor_bp']=='1000')
                credible=[r for r in bridges if r['site']==s['id'] and credible_bridge(r)]
                n=len({r['read'] for r in credible});unique=len({r['read'] for r in credible if min(int(r['unique_left_bp']),int(r['unique_right_bp']))>=1000})
                kind={'transition_midpoint':'transition midpoint','regional_trough':'support trough','gap_option':'scaffold gap'}[s['kind']]
                case_lines.append(f"| [{s['id']}]({assembly}/{e['event']}.{s['id']}.detail.png) ({kind}) | {s['point']:,} | {old['mapq20']} / {old['mapq30']} | {n} | {unique} |")
                actions.append(dict(assembly=assembly,id=s['id'],scaffold=e['scaffold'],selected='NO',action='BREAK_PROBABLE_MISJOIN',
                                    cut_bp='',assay_probe_bp=s['point'],reviewer='',reason='',coordinate_stage='pre_finishing',
                                    assessment_sha256=job['assessment_sha256'],decision_source='review',evidence_packet_id=key,
                                    localization_status='localized',auto_eligible='no'))
                summaries.append(dict(assembly=assembly,event=e['event'],site=s['id'],role=e['role'],position=s['point'],
                                      original_mapq20=old['mapq20'],original_mapq30=old['mapq30'],competitive_credible=n,unique_arm_bridges=unique))
            case_lines += ['', 'The transition midpoint is a measurement location, **not a proposed cut**. '
                           'Credible spans require MAPQ >=30, identity >=95%, at least 80% of the read aligned and no comparable same-haplotype alternative. Counts are descriptive, not a haplotype phasing test. '
                           'Other-haplotype similarity alone does not invalidate a bridge. Click a site ID for its detailed view.', '',
                           '### Can a read reach the unique chromosome anchors?','',
                           '| Site | Minimum span needed | Local read median / 95th percentile | Local reads long enough |',
                           '|---|---:|---:|---:|']
            for r in reach:
                if r['event']!=e['event']:continue
                case_lines.append(f"| {r['site']} | {r['minimum_unique_anchor_span_bp'] or 'No usable anchors'} | {r['local_read_median']} / {r['local_read_p95']} | {r['local_reads_long_enough'] or '0'} |")
            case_lines += ['', 'Read length feasibility is not an expected bridge count: start positions, coverage and phasing also matter.', '',
                           '### Chromosome-labelled islands beyond each site','',
                           '| Site | Chromosome on opposite side | Independent individuals with unique support | Largest unique block in a peer |',
                           '|---|---|---:|---:|']
            island_rows=table(dest/'labels'/'opposite_side_islands.tsv')
            for s in e['sites']:
                for chrom,side in ((e['left'],'right'),(e['right'],'left')):
                    rs=[r for r in island_rows if r['event']==e['event'] and r['site']==s['id'] and r['chrom']==chrom and r['side']==side]
                    independent={r['sample'] for r in rs if r['sample'] and r['sample']!=job['sample'] and r['eligible']=='True' and int(r['locus_unique_bp'])>0}
                    case_lines.append(f"| {s['id']} | {chrom}, {side} side | {len(independent)} | {max((int(r['largest_unique_block_bp']) for r in rs),default=0):,} bp |")
            case_lines += ['', 'These are observations, not a rule that a cut must follow the last alignment. '
                           'A long unique block needs a different explanation from a short multiply placed island. '
                           f'[All island lengths and per-peer evidence]({assembly}/labels/opposite_side_islands.tsv).','',
                           '### Aligned evidence','',f'![Evidence]({assembly}/{e["event"]}.evidence.png)','']
            fig,axes=plt.subplots(4,1,figsize=(15,12),sharex=True,gridspec_kw={'height_ratios':[2,1,1,1]},layout='constrained')
            peers=sorted({r['peer'] for r in tracks_data if r['event']==e['event']})
            ordered_chroms=sorted((e['left'],e['right']),key=lambda c:[int(v) if v.isdigit() else v for v in re.split(r'(\d+)',c)])
            colors=dict(zip(ordered_chroms,('#0072B2','#D55E00')))
            for i,peer in enumerate(peers):
                grouped=defaultdict(list)
                for r in tracks_data:
                    if r['event']!=e['event'] or r['peer']!=peer or r['mode']!='locus_unique' or not r['chrom']:continue
                    grouped[colors.get(r['chrom'],'#999999')].append((int(r['lo'])/1e6,(int(r['hi'])-int(r['lo']))/1e6))
                for color,intervals in grouped.items():axes[0].broken_barh(intervals,(i-.35,.7),facecolors=color)
            axes[0].set_yticks(range(len(peers)),peers);axes[0].set_title('Unique-locus chromosome assignments; blank = unresolved, not absent sequence',pad=45,fontsize=11)
            rs=[r for r in hifi if r['event']==e['event']]
            for tier,c in (('mapq20','#0072B2'),('mapq30','#D55E00')):
                axes[1].plot([int(r['position'])/1e6 for r in rs],[int(r[tier]) for r in rs],label=tier,color=c,lw=1)
            axes[1].set_ylabel('HiFi spanning\nmolecules');axes[1].legend()
            for lib in sorted({r['library'] for r in contacts}-{'POOLED'}):
                rs=[r for r in contacts if r['event']==e['event'] and r['library']==lib and r['flank_bp']=='300000']
                x=[int(r['point'])/1e6 for r in rs]
                y=[float(r['ratio']) if r['ratio'] else math.nan for r in rs]
                axes[2].plot(x,y,label=lib,lw=1)
                low=[i for i,r in enumerate(rs) if min(int(r['left_within']),int(r['right_within']))<100]
                axes[2].scatter([x[i] for i in low],[y[i] for i in low],marker='x',color='grey',s=12)
                axes[3].plot(x,[min(int(r['left_within']),int(r['right_within'])) for r in rs],label=lib,lw=1)
            axes[2].set_ylabel('Hi-C cross /\nâˆš(within L Ã— R)');axes[2].legend();axes[2].set_ylim(bottom=0)
            axes[3].set_ylabel('Weaker flank\nwithin-pairs');axes[3].set_yscale('symlog',linthresh=1);axes[3].set_ylim(bottom=0)
            for ax in axes:
                if ax is not axes[0]:
                    shade_masks(ax,masks[e['event']],colors,.08)
                ax.axvspan(e['transition_lo']/1e6,e['transition_hi']/1e6,color='#E69F00',alpha=.12)
                for s in e['sites']:ax.axvline(s['point']/1e6,color='#882255',ls='--',lw=.7,alpha=.7)
                ax.set_facecolor('white');ax.grid(alpha=.15)
            for i,s in enumerate(e['sites']):axes[0].text(s['point']/1e6,1.02+(i%2)*.07,s['id'],transform=axes[0].get_xaxis_transform(),ha='center',fontsize=7)
            axes[-1].set_xlabel('Assessment scaffold coordinate (Mb)');axes[-1].set_xlim(e['lo']/1e6,e['hi']/1e6)
            fig.suptitle(f"{assembly} / {e['event']}: {e['left']} â†’ {e['right']} Â· blue={ordered_chroms[0]}, orange={ordered_chroms[1]}")
            fig.savefig(dest/(e['event']+'.evidence.png'),dpi=140,facecolor='white');plt.close(fig)
            # Every measured site has a detailed view, even if it is outside the
            # chromosome-transition interval or has no recommended exact cut.
            for s in e['sites']:
                fig,axs=plt.subplots(3,1,figsize=(12,7),sharex=True,layout='constrained')
                for tier,c in (('mapq20','#0072B2'),('mapq30','#D55E00')):
                    rs=[r for r in hifi if r['event']==e['event']]
                    axs[0].plot([int(r['position'])/1e6 for r in rs],[int(r[tier]) for r in rs],color=c,label=tier,lw=1)
                axs[0].legend();axs[0].set_ylabel('HiFi molecules')
                for lib in sorted({r['library'] for r in contacts}-{'POOLED'}):
                    rs=[r for r in contacts if r['event']==e['event'] and r['library']==lib and r['flank_bp']=='100000']
                    axs[1].plot([int(r['point'])/1e6 for r in rs],[float(r['ratio']) if r['ratio'] else math.nan for r in rs],label=lib,lw=1)
                    axs[2].plot([int(r['point'])/1e6 for r in rs],[min(int(r['left_within']),int(r['right_within'])) for r in rs],label=lib,lw=1)
                axs[1].set_ylabel('Hi-C ratio\n100 kb flanks');axs[1].legend();axs[1].set_ylim(bottom=0);axs[2].set_ylabel('Weaker flank\nwithin-pairs')
                for ax in axs:
                    shade_masks(ax,masks[e['event']],colors,.12)
                    ax.axvline(s['point']/1e6,color='#882255',ls='--');ax.grid(alpha=.15);ax.set_facecolor('white')
                axs[-1].set_xlim(max(e['lo'],s['point']-100000)/1e6,min(e['hi'],s['point']+100000)/1e6)
                axs[-1].set_xlabel('Assessment scaffold coordinate (Mb)');fig.suptitle(f"{assembly} / {s['id']} Â· probe {s['point']:,} bp")
                fig.savefig(dest/f"{e['event']}.{s['id']}.detail.png",dpi=130,facecolor='white');plt.close(fig)
            case_lines += ['Grey crosses mark contact ratios with fewer than 100 within-flank pairs on either side; '
                           'this is an observability flag, not a significance threshold. '
                           'The fourth panel shows the actual weaker-flank count. Dashed lines mark every tested site.', '',
                           '### Unique chromosome-arm contacts','',
                           '| Site | Library | Unique sequence per arm | Cross pairs | Within left / right | Expected cross pairs from matched controls |',
                           '|---|---|---:|---:|---:|---:|']
            for r in anchors:
                if r['event']!=e['event'] or r['anchor_bp']!='100000' or r['library']=='POOLED':continue
                if r['status']!='measured':case_lines.append(f"| {r['site']} | {r['library']} | Insufficient unique sequence | â€” | â€” | â€” |")
                else:case_lines.append(f"| {r['site']} | {r['library']} | 100 kb | {r['cross']} | {r['left_within']} / {r['right_within']} | {round(float(r['expected_cross']),1) if r['expected_cross'] else 'Insufficient matched controls'} |")
            case_lines += ['', 'Selected unique bins can be discontinuous. Full coordinates and other anchor sizes are in '
                           f'[{assembly}/hic/unique_arm_contacts.tsv]({assembly}/hic/unique_arm_contacts.tsv). '
                           'Longer separation naturally reduces contacts; raw cross counts alone do not rank cuts.', '',
                           '### Local contact maps','']
            for record in maps:
                if record['event']!=e['event'] or record['library']=='POOLED':continue
                arr=np.load(dest/'hic'/record['file']);m=arr['matrix'];lo=record['lo'];res=record['resolution'];extent=[lo/1e6,(lo+len(m)*res)/1e6]*2
                fig,ax=plt.subplots(figsize=(7,6),layout='constrained')
                im=ax.imshow(np.log1p(m),origin='lower',extent=extent,cmap='magma',interpolation='nearest')
                for s in e['sites']:
                    ax.axvline(s['point']/1e6,color='cyan',lw=.5);ax.axhline(s['point']/1e6,color='cyan',lw=.5)
                ax.set_title(f"{e['event']} Â· {record['library']} Â· {res//1000} kb Â· raw counts")
                ax.set_xlabel('Assessment coordinate (Mb)');ax.set_ylabel('Assessment coordinate (Mb)');fig.colorbar(im,ax=ax,label='log(1 + pairs)')
                name=f"{e['event']}.{record['library']}.{res}.png";fig.savefig(dest/name,dpi=120,facecolor='white');plt.close(fig)
                if res==25000:case_lines += [f'![Contact map]({assembly}/{name})','']
                else:case_lines += [f'[{record["library"]}: {res//1000} kb contact map]({assembly}/{name})','']
            old_figure=Path(job['evidence']).parent/'plots'/(e['event']+'.png')
            if old_figure.exists():
                shutil.copyfile(old_figure,dest/(e['event']+'.original.png'))
                case_lines += ['<details>','<summary>Original overview, including depth and telomere evidence</summary>','',
                               f'![Original evidence overview]({assembly}/{e["event"]}.original.png)','','</details>','']
            case_lines += ['### Interpretation and limitations','',
                           'A break is favored by a chromosome-arm change plus loss of credible molecule continuity and '
                           'a contact deficit where both arms are observable. Retention is favored by multiple credible '
                           'bridging molecules with adequate anchors and no comparably good same-haplotype alternatives. '
                           'A trough entirely within one uniquely assigned arm may be a separate local problem.', '',
                           '**Decision remains for review. This diagnostic does not apply a new production scoring policy.**', '',
                           f'- [Molecule placements]({assembly}/hifi/bridge_placements.tsv)',
                           f'- [Chromosome evidence by 5 kb tile]({assembly}/labels/chromosome_tiles.tsv)',
                           f'- [Raw alignment locations and orientation]({assembly}/labels/alignment_hits.tsv)',
                           f'- [Contact scan]({assembly}/hic/contact_profile.tsv)',
                           f'- [Same-chromosome control measurements]({assembly}/hic/controls.tsv)',
                           f'- [Native graph context and limitations]({assembly}/labels/summary.json)','']
    lines += case_lines
    lines += ['## Manual cut worksheet','',
              '[manual-cut-worksheet.tsv](manual-cut-worksheet.tsv) contains **unselected** rows with blank cut coordinates. '
              'Assay probes are not preselected cuts. After review, enter the chosen exact coordinate, reviewer and reason. '
              'Use the production manual-review interface to validate and apply selected actions; this investigation never edits FASTA.', '',
              '## Scope','',
              'This packet audits alignment ambiguity using retained assembly alignments and competitive HiFi remapping. '
              'It does not constitute a complete repeat annotation or haplotype phasing analysis. '
              'Same-chromosome controls and overlapping scan windows are descriptive comparisons, not independent statistical replicates. '
              'Contact maps and underlying compressed count matrices are retained for further review.','']
    (a.out/'README.md').write_text('\n'.join(lines),encoding='utf-8')
    write(a.out/'sites.tsv',summaries,list(summaries[0]))
    write(a.out/'manual-cut-worksheet.tsv',actions,list(actions[0]))
    save(a.out/'provenance.json',manifest)


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    prep=sub.add_parser('prepare')
    for key in ('source','focus','final-dir','out'):prep.add_argument('--'+key,type=Path,required=True)
    branch=sub.add_parser('branch');branch.add_argument('--kind',choices=('labels','hifi','hic'),required=True)
    branch.add_argument('--job',type=Path,required=True);branch.add_argument('--out',type=Path,required=True)
    branch.add_argument('--threads',type=int,default=8);branch.add_argument('--masks',type=Path)
    rep=sub.add_parser('report');rep.add_argument('--manifest',type=Path,required=True)
    rep.add_argument('--packet',dest='packets',action='append',type=Path,required=True);rep.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.command=='prepare':prepare(a)
    elif a.command=='report':report(a)
    else:
        a.out.mkdir(exist_ok=False);job=json.loads(a.job.read_text())
        if a.kind=='hic':hic(job,a.out,a.threads,a.masks)
        elif a.kind=='hifi':hifi(job,a.out,a.threads,a.masks)
        else:labels(job,a.out,a.threads)
        save(a.out/'provenance.json',dict(assembly=job['assembly'],assessment_sha256=job['assessment_sha256'],
             coordinate_stage=job['coordinate_stage'],job_manifest_sha256=sha(a.job),
             script_sha256=sha(__file__),commands=COMMANDS,
             tools={name:subprocess.check_output([name,'--version'],text=True,errors='replace').splitlines()[0]
                    for name in ('samtools','minimap2')}))


if __name__=='__main__':main()
