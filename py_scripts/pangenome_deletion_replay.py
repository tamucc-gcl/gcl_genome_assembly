"""Read-only comparison of native and chromosome-isolated deletion filtering.

The guard changes filtering groups only. Source sequences and coordinates are never
split. Original query names are restored before the normal PAF filter/assignment.
"""
import argparse
from array import array
from collections import Counter, defaultdict
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET

from pangenome_source_audit import read_table, write_table, query_identity, union_bp, low_retention
from pangenome_assignment_replay import (CACTUS_COMMIT, LocalJob, initialize_cactus_tools,
    identity_aliases, select_queries, mapping_intervals, write_intervals)

POLICIES = ('native','chromosome_isolated','disabled')
ASSIGNMENTS = ('baseline','regional_1Mb')
GUARD_VERSION = 'chromosome-isolated-v1'


def lines(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    return opener(path,'rt',encoding='utf-8')


def sha256(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def check_tools():
    from cactus.shared.version import cactus_commit
    if cactus_commit != CACTUS_COMMIT or importlib.metadata.version('cactus') != '3.2.1':
        raise ValueError('Deletion replay requires the audited Cactus 3.2.1 image/commit')
    initialize_cactus_tools()
    result = dict(cactus_commit=cactus_commit,guard_version=GUARD_VERSION,
                  guard_source_sha256=sha256(__file__))
    for name in ('filter-paf-deletions','vg','rgfa-split'):
        binary = shutil.which(name)
        if not binary: raise ValueError('Required pinned tool missing: '+name)
        result[name+'_sha256'] = sha256(binary)
    return result


def node_owners(gfa, output):
    """Same chromosome membership as reference-boundary BFS, without per-node BFS.

    Rank-zero nodes are boundaries. Only non-reference nodes are unioned; reference
    nodes must NOT join adjacent chromosomes through a composite graph edge.
    Components bordering multiple chromosomes are explicit errors when used.
    """
    ids, names, chromosomes = {}, [], {}
    parent, rank, ref = array('i'), array('B'), array('i')
    with lines(gfa) as handle:
        for line in handle:
            if not line.startswith('S\t'): continue
            f = line.rstrip().split('\t'); name = f[1]
            if name in ids: raise ValueError('Duplicate graph node: '+name)
            idx = len(parent); ids[name] = idx; names.append(name)
            tags = {t[:2]:t[5:] for t in f[3:] if len(t)>5 and t[2]==':' and t[4]==':'}
            chrom = tags.get('SN') if tags.get('SR')=='0' else None
            if chrom: chromosomes.setdefault(chrom,len(chromosomes)+1)
            parent.append(idx); rank.append(0); ref.append(chromosomes[chrom] if chrom else 0)
    if not chromosomes: raise ValueError('SV rGFA lacks rank-zero reference chromosome tags')
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    def edges():
        with lines(gfa) as handle:
            for line in handle:
                if line.startswith('L\t'):
                    f = line.split('\t')
                    if f[1] not in ids or f[3] not in ids: raise ValueError('Unknown graph link node')
                    yield ids[f[1]],ids[f[3]]
    for a,b in edges():
        if ref[a] or ref[b]: continue
        a,b = root(a),root(b)
        if a==b: continue
        if rank[a]<rank[b]: a,b=b,a
        parent[b]=a
        if rank[a]==rank[b]: rank[a]+=1
    owner = array('i',[0])*len(parent)
    for a,b in edges():
        if bool(ref[a]) == bool(ref[b]): continue
        if ref[b]: a,b=b,a
        component = root(b); value = owner[component]
        owner[component] = ref[a] if value in (0,ref[a]) else -1
    labels = {value:name for name,value in chromosomes.items()}
    def rows():
        for i,name in enumerate(names):
            c = ref[i] or owner[root(i)]
            yield dict(node=name,reference_chromosome=labels.get(c,''),
                       status='anchored' if c>0 else 'conflicting' if c<0 else 'unanchored')
    write_table(output,rows(),['node','reference_chromosome','status'])
    return dict(nodes=len(parent),reference_chromosomes=len(chromosomes))


def select_paf(source, dest, keys, aliases):
    count = 0
    with lines(source) as src, open(dest,'w',encoding='utf-8') as out:
        for line in src:
            if query_identity(line.split('\t',1)[0],aliases) in keys:
                out.write(line); count+=1
    return count


def prepare(gfa,paf,unfiltered,filter_log,identities,sources,out,threads):
    import cactus
    from cactus.refmap.cactus_minigraph import minigraph_gfa_from_pansn
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    provenance=check_tools()
    config=Path(cactus.__file__).parent/'cactus_progressive_config.xml'
    mapping=ET.parse(config).getroot().find('graphmap')
    if mapping.get('collapse','none')!='none' or int(mapping.get('delFilter','-1'))<=0:
        raise ValueError('Replay requires default positive deletion filter and no collapse')
    ledger=read_table(identities);selected=select_queries(read_table(sources))
    if not selected: raise ValueError('No diagnostic sources selected')
    references=[r['graph_name'] for r in ledger if r['is_reference'].lower()=='true']
    if len(references)!=1: raise ValueError('Expected one audited reference')
    # Native sorting has no tie-breaker for equal query starts. Preserve the whole
    # original input and its order to reproduce that control exactly.
    pre_count=0
    with lines(unfiltered) as src,(out/'pre_deletion.paf').open('w',encoding='utf-8') as dst:
        for line in src: dst.write(line);pre_count+=1
    write_table(out/'selection.tsv',selected,list(selected[0]))
    shutil.copyfile(identities,out/'identities.tsv');shutil.copyfile(config,out/'effective_config.xml')
    shutil.copyfile(filter_log,out/'original_filter.log')
    job=LocalJob(out/'adapter_tmp')
    converted=minigraph_gfa_from_pansn(job,{r['graph_name'] for r in ledger},str(gfa),str(Path(gfa).resolve()))
    with gzip.open(converted,'rb') as src,(out/'cactus.sv.gfa').open('wb') as dst: shutil.copyfileobj(src,dst)
    ownership=node_owners(out/'cactus.sv.gfa',out/'node_owners.tsv.gz')
    cmd=['vg','convert','-r','0','-g',str(out/'cactus.sv.gfa'),'-p','-T',str(out/'graph.trans')]
    with (out/'graph.vg').open('wb') as dst: subprocess.run(cmd,stdout=dst,check=True)
    provenance.update(reference=references[0],selection_queries=len(selected),
        pre_deletion_records=pre_count,
        expected_native_sha256=sha256(paf),
        config_sha256=sha256(config),ownership=ownership,threads=threads,
        correction='native deletion filtering on separate (query,reference chromosome) lanes; restore original names',
        limitation='non-dominant chromosome lanes are also checked internally; this is a wrapper policy, not a rebuilt upstream executable')
    (out/'preparation.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    shutil.rmtree(job.root)


def lane_paf(source,owners_path,dest,ledger):
    owners={r['node']:(r['reference_chromosome'],r['status']) for r in read_table(owners_path)}
    lanes={};counts=Counter();matches=Counter()
    with lines(source) as src,open(dest,'w',encoding='utf-8') as out:
        for line in src:
            f=line.rstrip('\n').split('\t')
            if len(f)<12: raise ValueError('Malformed lane PAF')
            node=f[5].split('|',1)[1] if f[5].startswith('id=_MINIGRAPH_|') else f[5]
            if node not in owners or owners[node][1]!='anchored':
                raise ValueError('Cannot safely assign mapping node to one chromosome: '+node)
            key=f[0],owners[node][0]
            if key not in lanes: lanes[key]='gcl_lane_'+str(len(lanes)+1)
            f[0]=lanes[key];counts[key]+=1;matches[key]+=int(f[9])
            out.write('\t'.join(f)+'\n')
    write_table(ledger,(dict(lane=lane,query=key[0],reference_chromosome=key[1],
        records=counts[key],matching_bases=matches[key]) for key,lane in lanes.items()),
        ['lane','query','reference_chromosome','records','matching_bases'])
    return {lane:key[0] for key,lane in lanes.items()}


def restore_names(source,dest,names):
    with open(source,encoding='utf-8') as src,open(dest,'w',encoding='utf-8') as dst:
        for line in src:
            query,rest=line.split('\t',1)
            if query not in names: raise ValueError('Unknown filtering lane: '+query)
            dst.write(names[query]+'\t'+rest)


def deletion_command(prepared,source,config,threads):
    p=Path(prepared);m=ET.parse(config).getroot().find('graphmap')
    cmd=['filter-paf-deletions',str(p/'graph.vg'),str(p/'graph.trans'),str(source),
         '-d',m.get('delFilter'),'-v','-p','-t',str(threads)]
    for attr,flag in [('delFilterThreshold','-m'),('delFilterQuerySizeThreshold','-s')]:
        if m.get(attr) is not None: cmd += [flag,m.get(attr)]
    return cmd


def run_filter(prepared,policy,out,threads):
    from cactus.shared.configWrapper import ConfigWrapper
    from cactus.refmap.cactus_graphmap import filter_paf
    if policy not in POLICIES: raise ValueError('Unknown deletion policy')
    check_tools();p=Path(prepared).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    provenance=json.loads((p/'preparation.json').read_text(encoding='utf-8'))
    source=p/'pre_deletion.paf';post=out/'post_deletion.paf';cmd=[]
    if policy=='disabled': shutil.copyfile(source,post)
    else:
        if policy=='chromosome_isolated':
            names=lane_paf(source,p/'node_owners.tsv.gz',out/'lanes.paf',out/'lanes.tsv')
            source=out/'lanes.paf'
        cmd=deletion_command(p,source,p/'effective_config.xml',threads)
        target=out/'lane_filtered.paf' if policy=='chromosome_isolated' else post
        with target.open('w',encoding='utf-8') as dst,(out/'filter.log').open('w',encoding='utf-8') as err:
            subprocess.run(cmd,stdout=dst,stderr=err,check=True)
        if policy=='chromosome_isolated':
            restore_names(target,post,names);source.unlink();target.unlink()
    native_match=sha256(post)==provenance['expected_native_sha256'] if policy=='native' else None
    if policy=='native' and not native_match:
        raise ValueError('Native replay does not exactly reproduce the cached exported PAF; inspect logs before interpreting correction')
    selection=read_table(p/'selection.tsv');aliases=identity_aliases(read_table(p/'identities.tsv'))
    keys={(r['assembly'],r['sequence']) for r in selection}
    select_paf(post,out/'selected_post.paf',keys,aliases)
    post.unlink();(out/'selected_post.paf').rename(post)
    job=LocalJob(out/'adapter_tmp');config=ConfigWrapper(ET.parse(p/'effective_config.xml').getroot())
    result=filter_paf(job,str(post),config,reference=provenance['reference'])
    shutil.copyfile(result,out/'filtered.paf');shutil.rmtree(job.root)
    lengths={(r['assembly'],r['sequence']):int(r['source_bp']) for r in selection}
    stages={};summaries=[]
    for stage,path in [('pre_deletion',p/'pre_deletion.paf'),('post_deletion',post),('post_paf_filter',out/'filtered.paf')]:
        spans,counts=mapping_intervals(path,aliases,lengths);stages[stage]=spans
        for row in selection:
            key=row['assembly'],row['sequence']
            summaries.append(dict(policy=policy,assembly=key[0],sequence=key[1],stage=stage,
                source_bp=lengths[key],covered_bp=union_bp(spans[key]),records=counts[key]))
            if stage=='post_paf_filter': row['filtered_query_covered_bp']=union_bp(spans[key])
    write_table(out/'filter_coverage.tsv',summaries,list(summaries[0]))
    write_intervals(out/'coverage_intervals.tsv.gz',stages)
    write_table(out/'selection.tsv',selection,list(selection[0]))
    provenance.update(policy=policy,native_exact_match=native_match,command=cmd)
    (out/'preparation.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    for name in ('cactus.sv.gfa','effective_config.xml','identities.tsv'):
        (out/name).symlink_to((p/name).resolve())
    post.unlink()  # interval evidence survives; no duplicate large output


def report(filters,assignments,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    filter_rows=[];rows=[];native=False
    for directory in filters:
        p=Path(directory);v=json.loads((p/'preparation.json').read_text(encoding='utf-8'))
        if v['policy']=='native': native=v['native_exact_match'] is True
        filter_rows += read_table(p/'filter_coverage.tsv')
    if not native or {r['policy'] for r in filter_rows}!=set(POLICIES):
        raise ValueError('Incomplete or non-reproducing filter comparison')
    for directory in assignments:
        p=Path(directory);v=json.loads((p/'case.json').read_text(encoding='utf-8'))
        for row in read_table(p/'assignments.tsv'): rows.append(dict(policy=v['policy'],**row))
    keys=defaultdict(set)
    for r in rows: keys[r['assembly'],r['sequence']].add((r['policy'],r['case']))
    expected={(p,c) for p in POLICIES for c in ASSIGNMENTS}
    if not keys or any(v!=expected for v in keys.values()) or len(rows)!=len(keys)*len(expected):
        raise ValueError('Incomplete or duplicated assignment matrix')
    write_table(out/'assignment_comparison.tsv',rows,list(rows[0]))
    write_table(out/'filter_comparison.tsv',filter_rows,list(filter_rows[0]))
    index={(r['policy'],r['assembly'],r['sequence'],r['stage']):int(r['covered_bp']) for r in filter_rows}
    per={(r['policy'],r['assembly'],r['sequence'],r['case']):r for r in rows}
    selections=read_table(Path(filters[0])/'selection.tsv')
    text=['# Deletion-filter correction comparison','','Native filtering exactly reproduces the complete cached exported PAF byte for byte. '
        'These tests change intermediate mappings only; production assemblies and graphs are unchanged. Results require review before production use.','',
        'The chromosome-isolated guard runs the pinned deletion filter separately for each source/chromosome lane, '
        'then restores original query names. Each lane retains original query length and coordinates. '
        'It also checks non-dominant lanes internally, so this is a diagnostic wrapper policy, not an upstream binary patch.','',
        '| Source | Policy | Before deletion filter, Mb | After deletion filter, Mb | After PAF filter, Mb | Whole assignment, Mb | Regional assignment, Mb | Regional chromosomes |',
        '|---|---|---:|---:|---:|---:|---:|---|']
    for s in selections:
        if not (low_retention(s) or ('+' in s['sequence'] and int(s['source_bp'])>=10000000)): continue
        for policy in POLICIES:
            a,b=s['assembly'],s['sequence'];w=per[policy,a,b,'baseline'];r=per[policy,a,b,'regional_1Mb']
            numbers=[index[policy,a,b,stage] for stage in ('pre_deletion','post_deletion','post_paf_filter')]
            numbers += [int(w['surviving_query_covered_bp']),int(r['surviving_query_covered_bp'])]
            text.append(f"| {a} / {b} | {policy} | "+' | '.join(f'{v/1e6:.3f}' for v in numbers)+f" | {r['assigned_chromosomes'] or 'None'} |")
    deltas=[]
    text += ['', '## Cohort-wide mapping changes', '',
        '| Policy | Selected sources gaining coverage after PAF filtering | Sources losing coverage | Net change, Mb |',
        '|---|---:|---:|---:|']
    for policy in ('chromosome_isolated','disabled'):
        values=[]
        for s in selections:
            a,b=s['assembly'],s['sequence']
            old=index['native',a,b,'post_paf_filter'];new=index[policy,a,b,'post_paf_filter']
            values.append(new-old)
            deltas.append(dict(policy=policy,assembly=a,sequence=b,native_bp=old,alternative_bp=new,delta_bp=new-old))
        text.append(f'| {policy} | {sum(v>0 for v in values)} | {sum(v<0 for v in values)} | {sum(values)/1e6:.3f} |')
    write_table(out/'mapping_changes.tsv',deltas,['policy','assembly','sequence','native_bp','alternative_bp','delta_bp'])
    controls=[s for s in selections if s['chromosome_scale']=='yes' and '+' not in s['sequence'] and not low_retention(s)
              and s['status']=='represented_in_clip']
    regressions=[]
    for s in controls:
        a,b=s['assembly'],s['sequence']
        for case in ASSIGNMENTS:
            old=int(per['native',a,b,case]['surviving_query_covered_bp'])
            new=int(per['chromosome_isolated',a,b,case]['surviving_query_covered_bp'])
            if new<old: regressions.append(dict(assembly=a,sequence=b,case=case,native_bp=old,guarded_bp=new))
    write_table(out/'control_regressions.tsv',regressions,['assembly','sequence','case','native_bp','guarded_bp'])
    text += ['',f'{len(controls)} ordinary chromosome controls tested; {len(regressions)} control/case combinations have less surviving mapped coverage under the guard.', '',
        '[All filtering results](filter_comparison.tsv) · [All assignment results](assignment_comparison.tsv) · '
        '[Every source gain/decrease](mapping_changes.tsv) · [Control decreases](control_regressions.tsv)', '',
        '[Native provenance and command](review_native/preparation.json) · '
        '[Guard provenance and command](review_chromosome_isolated/preparation.json) · '
        '[Guard lane ledger](review_chromosome_isolated/lanes.tsv) · '
        '[Native filter log](review_native/filter.log) · '
        '[Guard filter log](review_chromosome_isolated/filter.log)', '',
        'Per-policy folders contain commands, tool hashes, filtering logs, lane ledgers and merged coverage intervals. '
        'Per-assignment folders contain decisions, region coordinates and surviving mapping intervals. '
        'Query coverage is an interval union, not proof of correct homology or a final graph path. '
        'Recovered arms must have credible chromosome-specific assignment; unexplained control decreases need review. '
        'Do not enable this policy in production or rebuild Cactus until these comparisons are accepted.','']
    (out/'deletion_filter_report.md').write_text('\n'.join(text),encoding='utf-8')


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='mode',required=True)
    prep=s.add_parser('prepare')
    for k in ('gfa','paf','unfiltered','filter_log','identities','sources','out'): prep.add_argument('--'+k.replace('_','-'),required=True)
    prep.add_argument('--threads',type=int,default=1)
    run=s.add_parser('filter')
    for k in ('prepared','policy','out'): run.add_argument('--'+k,required=True)
    run.add_argument('--threads',type=int,default=1)
    summary=s.add_parser('report');summary.add_argument('--filters',nargs='+',required=True)
    summary.add_argument('--assignments',nargs='+',required=True);summary.add_argument('--out',required=True)
    a=vars(p.parse_args());mode=a.pop('mode')
    {'prepare':prepare,'filter':run_filter,'report':report}[mode](**a)


if __name__=='__main__': main()
