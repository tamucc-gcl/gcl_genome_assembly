"""Bounded existing-option comparison for cached regional chromosome assignments."""
import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET

from pangenome_assignment_replay import split_command, summarize_split, identity_aliases, fragment, merged, interval_difference
from pangenome_regional_patch import sha, coordinate_audit, audit_paf_membership, PATCH
from pangenome_regional_trace import TRACE_VERSION
from pangenome_source_audit import read_table, write_table, union_bp, low_retention

CHUNKS=(1000000,100000,25000)
VARIANTS=('native','fixed')
POLICIES=('native','chromosome_isolated','disabled')


def chunk_command(gfa,paf,config,reference,contigs,chunk,prefix,log):
    if chunk not in CHUNKS:raise ValueError('Unsupported chunk comparison setting')
    cmd=split_command(gfa,paf,config,reference,contigs,'regional_1Mb',prefix,log)
    if cmd.count('-u')!=1:raise ValueError('Expected exactly one regional chunk option')
    cmd[cmd.index('-u')+1]=str(chunk)
    return cmd


def run(prepared,tools,fixed_binary,fixed_info,variant,chunk,contigs,out):
    if variant not in VARIANTS or chunk not in CHUNKS:raise ValueError('Unsupported chunk comparison case')
    p=Path(prepared);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    info=json.loads((Path(tools)/'trace_build.json').read_text());fixed=json.loads(Path(fixed_info).read_text())
    if info['version']!=TRACE_VERSION or fixed['patch']!=PATCH:raise ValueError('Unexpected comparison tool version')
    binary=Path(tools)/'native.bin' if variant=='native' else Path(fixed_binary)
    expected=info['binaries']['native'] if variant=='native' else fixed['patched_binary_sha256']
    if sha(binary)!=expected:raise ValueError('Comparison executable checksum mismatch')
    provenance=json.loads((p/'preparation.json').read_text());selection=read_table(p/'selection.tsv')
    ambiguous=ET.parse(p/'effective_config.xml').getroot().find('graphmap_split').get('ambiguousName','_AMBIGUOUS_')
    temp=out/'split_work';temp.mkdir()
    cmd=chunk_command(p/'cactus.sv.gfa',p/'filtered.paf',p/'effective_config.xml',provenance['reference'],
        contigs,chunk,temp/'split_',out/'assignment_log.txt');cmd[0]=str(binary.resolve())
    (out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    with (out/'tool_stderr.txt').open('w') as err:subprocess.run(cmd,check=True,stdout=err,stderr=err)
    case=f'chunk_{chunk}_{variant}'
    assignments,regions=summarize_split(temp,p/'selection.tsv',p/'identities.tsv',case,ambiguous,out/'surviving_intervals.tsv.gz')
    audits=coordinate_audit(regions,selection,ambiguous);audit_paf_membership(temp,regions,selection,p/'identities.tsv')
    lengths={(r['assembly'],r['sequence']):int(r['source_bp']) for r in selection}
    aliases=identity_aliases(read_table(p/'identities.tsv'));spans=defaultdict(list)
    for file in temp.glob('split_*.paf'):
        chrom=file.name[len('split_'):-len('.paf')]
        if chrom==ambiguous:continue
        with file.open() as f:
            for line in f:
                fields=line.rstrip().split('\t');parsed=fragment(fields[0],aliases,lengths)
                if parsed is None:continue
                key,offset,end=parsed;spans[(*key,chrom)].append((offset+int(fields[2]),offset+int(fields[3])))
    mapped=[dict(assembly=k[0],sequence=k[1],chromosome=k[2],mapped_union_bp=union_bp(v)) for k,v in sorted(spans.items())]
    write_table(out/'chromosome_mapping.tsv',mapped,['assembly','sequence','chromosome','mapped_union_bp'])
    write_table(out/'assignments.tsv',assignments,list(assignments[0]))
    write_table(out/'regions.tsv',regions,['case','assembly','sequence','chromosome','start_0based','end_exclusive'])
    write_table(out/'coordinate_audit.tsv',audits,list(audits[0]))
    shutil.copyfile(p/'selection.tsv',out/'selection.tsv')
    (out/'case.json').write_text(json.dumps(dict(policy=provenance['policy'],variant=variant,chunk_bp=chunk,
        ambiguous_name=ambiguous,coordinate_integrity='PASS',binary_sha256=expected))+'\n')
    shutil.rmtree(temp)


def key(r):return r['assembly'],r['sequence']


def region_index(rows,ambiguous):
    spans=defaultdict(list)
    for r in rows:
        if r['chromosome']!=ambiguous:spans[(*key(r),r['chromosome'])].append((int(r['start_0based']),int(r['end_exclusive'])))
    return {k:merged(v) for k,v in spans.items()}


def audit_rows(rows,selection):
    sources={key(r):int(r['source_bp']) for r in selection};observed={key(r):r for r in rows}
    if len(observed)!=len(rows) or observed.keys()!=sources.keys():raise ValueError('Incomplete coordinate audit sources')
    for k,r in observed.items():
        if (r['status']!='PASS' or int(r['duplicate_or_overlapping_bp'])!=0 or
            int(r['assigned_union_bp'])!=int(r['assigned_sum_bp']) or
            int(r['assigned_union_bp'])+int(r['ambiguous_bp'])+int(r['unreported_bp'])!=sources[k]):
            raise ValueError('Invalid coordinate accounting')


def load_cases(baselines,comparisons):
    cases={};selection=None;source_ids=None;ambiguous=None
    for directory in comparisons:
        p=Path(directory);meta=json.loads((p/'case.json').read_text());a=meta['ambiguous_name']
        if ambiguous is not None and a!=ambiguous:raise ValueError('Inconsistent ambiguous names')
        ambiguous=a
    if ambiguous is None:raise ValueError('Missing chunk comparisons')
    for kind,dirs in (('baseline',baselines),('comparison',comparisons)):
        for directory in dirs:
            p=Path(directory);sources=read_table(p/'selection.tsv')
            ids={key(r):int(r['source_bp']) for r in sources}
            if len(ids)!=len(sources):raise ValueError('Duplicate comparison source')
            if selection is None:selection=sources;source_ids=ids
            elif ids!=source_ids:raise ValueError('Inconsistent comparison source selection')
            meta=json.loads((p/('trace_case.json' if kind=='baseline' else 'case.json')).read_text())
            policy=meta['policy']
            if kind=='baseline':
                eq=read_table(p/'observation_equivalence.tsv')
                if len(eq)!=2 or {r['variant'] for r in eq}!=set(VARIANTS) or any(
                    r['byte_identical']!='yes' or r['coordinate_integrity']!='PASS' for r in eq):
                    raise ValueError('Unverified cached baseline')
            for variant in (VARIANTS if kind=='baseline' else (meta['variant'],)):
                chunk=1000000 if kind=='baseline' else int(meta['chunk_bp']);k=policy,variant,chunk
                if policy not in POLICIES or variant not in VARIANTS or chunk not in CHUNKS or k in cases:
                    raise ValueError('Invalid/duplicate chunk matrix case')
                prefix=variant+'_' if kind=='baseline' else ''
                assignments=read_table(p/(prefix+'assignments.tsv'));regions=read_table(p/(prefix+'regions.tsv'))
                audit_rows(read_table(p/(prefix+'coordinate_audit.tsv')),sources)
                coordinate_audit(regions,sources,ambiguous)
                assigned={key(r):r for r in assignments}
                if len(assigned)!=len(assignments) or assigned.keys()!=source_ids.keys():raise ValueError('Incomplete assignment sources')
                cases[k]=dict(assignments=assigned,regions=region_index(regions,ambiguous),
                    mapped={} if kind=='baseline' else {(*key(r),r['chromosome']):int(r['mapped_union_bp'])
                        for r in read_table(p/'chromosome_mapping.tsv')})
    if set(cases)!={(p,v,c) for p in POLICIES for v in VARIANTS for c in CHUNKS}:
        raise ValueError('Incomplete 18-case chunk comparison matrix')
    return cases,selection


def report(baselines,comparisons,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);cases,selection=load_cases(baselines,comparisons)
    controls={key(r) for r in selection if r['chromosome_scale']=='yes' and r['status']=='represented_in_clip'
              and '+' not in r['sequence'] and not low_retention(r)}
    sources={key(r):r for r in selection};rows=[];chrom_rows=[];summaries=[]
    fields=['policy','variant','chunk_bp','assembly','sequence','chromosome','change','start_0based','end_exclusive']
    with gzip.open(out/'chunk_region_changes.tsv.gz','wt',encoding='utf-8',newline='') as handle:
        import csv
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for (policy,variant,chunk),case in sorted(cases.items()):
            baseline=cases[policy,variant,1000000];changed=[];new_chroms=[]
            for source,r in sorted(case['assignments'].items()):
                before=baseline['assignments'][source]
                old_chroms={k[2] for k in baseline['regions'] if k[:2]==source}
                new_chromosomes={k[2] for k in case['regions'] if k[:2]==source}
                delta=int(r['surviving_query_covered_bp'])-int(before['surviving_query_covered_bp'])
                new=new_chromosomes-old_chroms;lost=old_chroms-new_chromosomes
                row=dict(policy=policy,variant=variant,chunk_bp=chunk,assembly=source[0],sequence=source[1],
                    ordinary_control='yes' if source in controls else 'no',baseline_mapped_bp=before['surviving_query_covered_bp'],
                    mapped_bp=r['surviving_query_covered_bp'],mapped_delta_bp=delta,baseline_chromosomes=','.join(sorted(old_chroms)),
                    chromosomes=','.join(sorted(new_chromosomes)),new_chromosomes=','.join(sorted(new)),lost_chromosomes=','.join(sorted(lost)),
                    assigned_span_bp=r['assigned_span_bp'],ambiguous_bp=r['ambiguous_span_bp'])
                rows.append(row);changed.append(row)
                for chrom in sorted(old_chroms|new_chromosomes):
                    k=(*source,chrom);old=baseline['regions'].get(k,[]);spans=case['regions'].get(k,[])
                    added=interval_difference(spans,old);removed=interval_difference(old,spans)
                    for label,intervals in (('added',added),('removed',removed)):
                        for lo,hi in intervals:writer.writerow(dict(policy=policy,variant=variant,chunk_bp=chunk,assembly=source[0],
                            sequence=source[1],chromosome=chrom,change=label,start_0based=lo,end_exclusive=hi))
                    if not added and not removed:continue
                    largest=max((hi-lo for lo,hi in spans),default=0);assigned_bp=union_bp(spans)
                    cr=dict(policy=policy,variant=variant,chunk_bp=chunk,assembly=source[0],sequence=source[1],chromosome=chrom,
                        newly_assigned='yes' if chrom in new else 'no',baseline_assigned_bp=union_bp(old),assigned_bp=assigned_bp,
                        added_bp=union_bp(added),removed_bp=union_bp(removed),region_count=len(spans),largest_region_bp=largest,
                        mapped_union_bp=case['mapped'].get(k,0),mapped_fraction=case['mapped'].get(k,0)/assigned_bp if assigned_bp else 0,
                        ordinary_control='yes' if source in controls else 'no',
                        review_flag='new chromosome on ordinary control' if chrom in new and source in controls else
                            'new chromosome represented only by regions below 1 Mb' if chrom in new and largest<1000000 else '')
                    chrom_rows.append(cr)
                    if chrom in new:new_chroms.append(cr)
            summaries.append(dict(policy=policy,variant=variant,chunk_bp=chunk,sources=len(changed),
                mapped_bp=sum(int(r['mapped_bp']) for r in changed),mapped_delta_bp=sum(int(r['mapped_delta_bp']) for r in changed),
                sources_gaining=sum(int(r['mapped_delta_bp'])>0 for r in changed),sources_losing=sum(int(r['mapped_delta_bp'])<0 for r in changed),
                controls_losing=sum(r['ordinary_control']=='yes' and int(r['mapped_delta_bp'])<0 for r in changed),
                controls_with_new_chromosome=sum(r['ordinary_control']=='yes' and bool(r['new_chromosomes']) for r in changed),
                new_chromosome_assignments=len(new_chroms),new_chromosomes_below_1Mb=sum(int(r['largest_region_bp'])<1000000 for r in new_chroms)))
    write_table(out/'chunk_source_comparison.tsv',rows,list(rows[0]))
    write_table(out/'chunk_matrix_summary.tsv',summaries,list(summaries[0]))
    chrom_fields=['policy','variant','chunk_bp','assembly','sequence','chromosome','newly_assigned','baseline_assigned_bp',
        'assigned_bp','added_bp','removed_bp','region_count','largest_region_bp','mapped_union_bp','mapped_fraction','ordinary_control','review_flag']
    write_table(out/'chunk_chromosome_changes.tsv',chrom_rows,chrom_fields)
    flagged=[r for r in chrom_rows if r['review_flag'] or int(r['removed_bp'])>0]
    write_table(out/'chunk_review_flags.tsv',flagged,chrom_fields)
    lines=['# Regional chunk-setting comparison','',
        'Only the existing `-u` option changes: 1 Mb, 100 kb and 25 kb. The 1 Mb baselines come from the verified cached gate run. '
        'Both original and corrected grouping are compared under all three filtering policies. Production graph defaults and assembly cuts are unchanged.','',
        '| Filtering | Grouping | Chunk, kb | Mapped change, Mb | Sources gaining / losing | Controls losing / gaining chromosomes | New chromosome assignments / only small regions |',
        '|---|---|---:|---:|---|---|---|']
    for r in summaries:
        lines.append(f"| {r['policy']} | {r['variant']} | {r['chunk_bp']/1000:g} | {r['mapped_delta_bp']/1e6:+.3f} | "
            f"{r['sources_gaining']} / {r['sources_losing']} | {r['controls_losing']} / {r['controls_with_new_chromosome']} | "
            f"{r['new_chromosome_assignments']} / {r['new_chromosomes_below_1Mb']} |")
    lines+=['','## Large composites and low-retention sources','',
        '| Source | Filtering | Grouping | Chunk, kb | Mapped, Mb | Change, Mb | Assigned chromosomes | Ambiguous, Mb |',
        '|---|---|---|---:|---:|---:|---|---:|']
    for r in rows:
        source=sources[key(r)]
        if not (low_retention(source) or ('+' in source['sequence'] and int(source['source_bp'])>=10000000)):continue
        lines.append(f"| {r['assembly']} / {r['sequence']} | {r['policy']} | {r['variant']} | {r['chunk_bp']/1000:g} | "
            f"{int(r['mapped_bp'])/1e6:.3f} | {r['mapped_delta_bp']/1e6:+.3f} | {r['chromosomes'] or 'None'} | {int(r['ambiguous_bp'])/1e6:.3f} |")
    lines+=['','## Side effects and ownership','',
        f'{len(selection)} selected sources in each of 18 cases; {len(controls)} ordinary chromosome controls. '
        'Every case must pass source-bound and non-overlap checks, and new runs also verify each output mapping belongs to its assigned fragment. '
        'Reported mapped coverage uses source-coordinate unions.','',
        'Review new chromosome assignments on controls, mapping losses, and newly assigned chromosomes represented only by regions below 1 Mb. '
        'These flags identify side effects for inspection; they do not establish biological errors. The complete selected cohort includes fragmented assemblies.','',
        '[All sources](chunk_source_comparison.tsv) · [Chromosome/region changes](chunk_chromosome_changes.tsv) · '
        '[Review flags](chunk_review_flags.tsv) · [Exact added/removed intervals](chunk_region_changes.tsv.gz)','',
        'This is a bounded parameter experiment. Recovery of the motivating examples alone is insufficient to select a production default.','']
    (out/'regional_chunk_report.md').write_text('\n'.join(lines),encoding='utf-8')


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='mode',required=True)
    r=s.add_parser('run')
    for k in ('prepared','tools','fixed_binary','fixed_info','variant','out'):r.add_argument('--'+k.replace('_','-'),required=True)
    r.add_argument('--chunk',type=int,choices=CHUNKS,required=True);r.add_argument('--contigs',nargs='+',required=True)
    r=s.add_parser('report');r.add_argument('--baselines',nargs='+',required=True);r.add_argument('--comparisons',nargs='+',required=True);r.add_argument('--out',required=True)
    args=vars(p.parse_args());mode=args.pop('mode');{'run':run,'report':report}[mode](**args)


if __name__=='__main__':main()
